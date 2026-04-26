import config

import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

from data_utils import read_10x_mtx, write_10x_mtx

adata = read_10x_mtx(config.hepatocytes_count_dir)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'])
umap_df['cell'] = adata.obs_names

# K-means clustering on UMAP coordinates
kmeans = KMeans(n_clusters=config.N_KMEANS_CLUSTERS, random_state=42, n_init='auto')
umap_df['kmeans_cluster'] = kmeans.fit_predict(umap_df[['UMAP_1', 'UMAP_2']])

# Plot
plt.rcParams['savefig.bbox'] = 'tight'
plt.figure(figsize=(10, 8))
sns.set_style("white")
palette = sns.color_palette('tab10', n_colors=config.N_KMEANS_CLUSTERS)
ax = sns.scatterplot(
    data=umap_df,
    x='UMAP_1',
    y='UMAP_2',
    hue='kmeans_cluster',
    palette=palette,
    s=1,
    edgecolor=None,
    alpha=0.7,
    legend=False
)
for cid, group in umap_df.groupby('kmeans_cluster'):
    ax.text(group['UMAP_1'].median(), group['UMAP_2'].median(), str(cid),
            fontsize=12, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.6, ec='none'))

plt.title(f"Hepatocytes: K-means clusters (k={config.N_KMEANS_CLUSTERS})", fontsize=15)
plt.axis('off')

os.makedirs(config.OUTPUT_DIR, exist_ok=True)
output_file = os.path.join(config.OUTPUT_DIR, "hepatocyte_umap")
plt.savefig(output_file, dpi=300)

umap_csv = os.path.join(config.OUTPUT_DIR, "umap_hepatocytes_kmeans.csv")
umap_df.to_csv(umap_csv, index=False)
print(f"Saved UMAP CSV -> {umap_csv}")

# Write a matrix for each cluster
adata.obs['kmeans_cluster'] = umap_df['kmeans_cluster'].values
for cid in range(config.N_KMEANS_CLUSTERS):
    cluster_adata = adata[adata.obs['kmeans_cluster'] == cid]
    out_dir = os.path.join(config.hepatocytes_clusters_dir, f"cluster_{cid}")
    write_10x_mtx(cluster_adata, out_dir)
