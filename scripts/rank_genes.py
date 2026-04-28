import config

import os
import pandas as pd
import scanpy as sc

from data_utils import read_10x_mtx

# Load hepatocyte counts and attach cluster labels
adata = read_10x_mtx(config.hepatocytes_count_dir)

umap_df = pd.read_csv(os.path.join(config.OUTPUT_DIR, "umap_hepatocytes_kmeans.csv"))
umap_df = umap_df.set_index("cell")
adata.obs["cluster"] = umap_df.loc[adata.obs_names, "kmeans_cluster"].astype(str).values

# Normalize and log-transform before ranking
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Wilcoxon rank-sum test: each cluster vs. the rest
sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", key_added="wilcoxon")

# Build a tidy CSV: one row per (cluster, gene)
frames = []
for cluster in adata.obs["cluster"].unique():
    result = sc.get.rank_genes_groups_df(adata, group=cluster, key="wilcoxon")
    result.insert(0, "cluster", cluster)
    frames.append(result)

df = pd.concat(frames).reset_index(drop=True)
df = df.rename(columns={"names": "gene", "scores": "score", "logfoldchanges": "log2fc",
                         "pvals": "pval", "pvals_adj": "pval_adj"})
df = df[["cluster", "gene", "score", "log2fc", "pval", "pval_adj"]]

os.makedirs(config.OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(config.OUTPUT_DIR, "cluster_marker_genes.csv")
df.to_csv(out_path, index=False)
print(f"Saved {len(df)} rows -> {out_path}")
print(df.groupby("cluster").head(3)[["cluster", "gene", "score", "log2fc", "pval_adj"]].to_string(index=False))
