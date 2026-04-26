import config

import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from data_utils import read_10x_mtx

# Read the highly variable count matrix
adata = read_10x_mtx(config.highly_variable_count_dir)

# Attach cell type annotations from the annotation file
annot_df = pd.read_csv(config.annotated_liver_path, index_col='cell')
adata.obs['annot'] = annot_df.loc[adata.obs_names, 'annot'].values

# Normalize and compute UMAP
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot
plt.rcParams['savefig.bbox'] = 'tight'
umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'])
umap_df['annot'] = adata.obs['annot'].values

plt.figure(figsize=(10, 8))
sns.set_style("white")
ax = sns.scatterplot(
    data=umap_df,
    x='UMAP_1',
    y='UMAP_2',
    hue='annot',
    palette='tab20',
    s=1,
    edgecolor=None,
    alpha=0.7
)
plt.title("Liver Atlas: Highly Variable Genes", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Cell Type', markerscale=10)
plt.axis('off')

output_file = os.path.join(config.OUTPUT_DIR, "highly_variable_umap")
plt.savefig(output_file, dpi=300)
plt.show()
