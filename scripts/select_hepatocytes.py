import config

import pandas as pd
from data_utils import read_10x_mtx, write_10x_mtx

adata = read_10x_mtx(config.highly_variable_count_dir)

# Attach annotations and filter to hepatocytes
annot_df = pd.read_csv(config.annotated_liver_path, index_col='cell')
adata.obs['annot'] = annot_df.loc[adata.obs_names, 'annot'].values

adata = adata[adata.obs['annot'] == 'Hepatocytes']
print(f"Kept:    {adata.n_obs} hepatocytes")

write_10x_mtx(adata, config.hepatocytes_count_dir)
