import config

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#read the file as a dataframe
df = pd.read_csv(config.annotated_liver_path)

# 2. Set the style
plt.rcParams['savefig.bbox'] = 'tight'
plt.figure(figsize=(10, 8))
sns.set_style("white")

# 3. Draw the UMAP
# We color it by 'annot' (the cell type)
ax = sns.scatterplot(
    data=df, 
    x='UMAP_1', 
    y='UMAP_2', 
    hue='annot', 
    palette='tab20', 
    s=1,          
    edgecolor=None, 
    alpha=0.7
)

# 4. Clean up the plot
plt.title("Liver Atlas: Human CD45-negative Cells", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Cell Type', markerscale=10)
plt.axis('off')

# 5. write the plot as a png
output_file = os.path.join(config.OUTPUT_DIR, "liver_atlas_umap")
plt.savefig(output_file, dpi=300)

