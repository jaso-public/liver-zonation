import config
import scanpy as sc

from data_utils import read_10x_mtx, write_10x_mtx

adata = read_10x_mtx(config.selected_cells_count_dir)

# Select top N_HVG highly variable genes using raw counts (pearson residuals)
sc.experimental.pp.highly_variable_genes(adata, n_top_genes=config.N_HVG, flavor='pearson_residuals')
adata = adata[:, adata.var.highly_variable]
print(f"Kept:    {adata.n_vars} highly variable genes")

write_10x_mtx(adata, config.highly_variable_count_dir)
