import config

import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

from data_utils import read_10x_mtx

# Load hepatocyte counts and attach cluster labels
adata = read_10x_mtx(config.hepatocytes_count_dir)

umap_df = pd.read_csv(os.path.join(config.OUTPUT_DIR, "umap_hepatocytes_kmeans.csv"))
umap_df = umap_df.set_index("cell")
adata.obs["cluster"] = umap_df.loc[adata.obs_names, "kmeans_cluster"].astype(str).values

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", key_added="wilcoxon")

os.makedirs(config.OUTPUT_DIR, exist_ok=True)

# ── 1. Dot plot: fraction expressing × mean expression per cluster
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=config.N_TOP_MARKERS,
    key="wilcoxon",
    standard_scale="var",
    show=False,
)
out = os.path.join(config.OUTPUT_DIR, "marker_genes_dotplot.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved -> {out}")

# ── 2. Log2FC bar charts: top N marker genes per cluster
markers = pd.read_csv(os.path.join(config.OUTPUT_DIR, "cluster_marker_genes.csv"))
top = (markers[markers["pval_adj"] < 0.05]
       .sort_values("score", ascending=False)
       .groupby("cluster")
       .head(config.N_TOP_MARKERS))

clusters = sorted(top["cluster"].unique())
palette = sns.color_palette("tab10", n_colors=len(clusters))

fig, axes = plt.subplots(1, len(clusters), figsize=(4 * len(clusters), 5), sharey=False)
if len(clusters) == 1:
    axes = [axes]

for ax, cid, color in zip(axes, clusters, palette):
    grp = top[top["cluster"] == cid].sort_values("log2fc")
    ax.barh(grp["gene"], grp["log2fc"], color=color)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_title(f"Cluster {cid}", fontsize=11)
    ax.set_xlabel("log2FC")
    if ax is axes[0]:
        ax.set_ylabel("Gene")

fig.suptitle(f"Top {config.N_TOP_MARKERS} marker genes per cluster (Wilcoxon, FDR < 0.05)", y=1.02)
plt.tight_layout()
out = os.path.join(config.OUTPUT_DIR, "marker_genes_log2fc.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved -> {out}")
