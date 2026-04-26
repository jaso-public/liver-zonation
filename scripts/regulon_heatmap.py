import config

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

# Load and concatenate AUC matrices from all clusters
frames = []
for cluster_id in range(config.N_KMEANS_CLUSTERS):
    auc_path = os.path.join(config.hepatocytes_scenic_dir, f"cluster_{cluster_id}", "auc_matrix.csv")
    if not os.path.exists(auc_path):
        print(f"Cluster {cluster_id}: AUC matrix not found, skipping")
        continue
    auc = pd.read_csv(auc_path, index_col=0)
    auc['cluster'] = cluster_id
    frames.append(auc)

df = pd.concat(frames).fillna(0)
cluster_labels = df.pop('cluster')

# Strip pySCENIC "(+)" suffix from column names
df.columns = [c.split('(')[0].strip() for c in df.columns]

# Keep only top N_TOP_TFS regulons by mean activity
top_regulons = df.mean().nlargest(config.N_TOP_TFS).index
df = df[top_regulons]

# Within each cluster, sort cells by hierarchical clustering
ordered_idx = []
for cid in range(config.N_KMEANS_CLUSTERS):
    mask = cluster_labels == cid
    sub = df[mask]
    if len(sub) < 2:
        ordered_idx.extend(sub.index.tolist())
        continue
    Z = linkage(sub.values, method='ward')
    order = leaves_list(Z)
    ordered_idx.extend(sub.iloc[order].index.tolist())

df_ordered = df.loc[ordered_idx]
cluster_ordered = cluster_labels.loc[ordered_idx]

# Row color bar by cluster
palette = sns.color_palette('tab10', n_colors=config.N_KMEANS_CLUSTERS)
row_colors = cluster_ordered.map(lambda c: palette[c])

# Plot
os.makedirs(config.OUTPUT_DIR, exist_ok=True)
g = sns.clustermap(
    df_ordered,
    row_cluster=False,
    col_cluster=True,
    row_colors=row_colors,
    cmap='viridis',
    xticklabels=True,
    yticklabels=False,
    figsize=(12, 10),
    cbar_kws={'label': 'AUC'},
)
g.ax_heatmap.set_xlabel("Regulon")
g.ax_heatmap.set_ylabel("Cells")
g.figure.suptitle(f"Regulon activity by cluster (top {config.N_TOP_TFS} TFs)", y=1.02, fontsize=13)

handles = [plt.Rectangle((0, 0), 1, 1, color=palette[c], label=f"Cluster {c}")
           for c in range(config.N_KMEANS_CLUSTERS)]
g.ax_heatmap.legend(handles=handles, bbox_to_anchor=(1.25, 1), loc='upper left', title='Cluster')

out_path = os.path.join(config.OUTPUT_DIR, "regulon_heatmap")
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.show()
print(f"Saved -> {out_path}")
