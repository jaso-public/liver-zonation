import config

import pandas as pd
from data_utils import read_10x_mtx, write_10x_mtx

#read the file as a dataframe
df = pd.read_csv(config.annotated_liver_path)

# Now we read the count table as an AnnData object.
adata = read_10x_mtx(config.count_table_human_dir)

# Now we subset the adata to only keep the cells that are in the annotation file.
annot_set = set(df['cell'])
adata = adata[adata.obs_names.isin(annot_set)]
print(f"Kept:    {adata.n_obs} annotated cells ({len(annot_set)} in annotation file)")

# write the subsetted adata as a new count table in the output directory
write_10x_mtx(adata, config.selected_cells_count_dir)
