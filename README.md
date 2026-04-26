# Hepatocyte Regulon Analysis
Liver atlas single-cell pipeline: HVG selection → hepatocyte clustering → pySCENIC regulon inference → visualization.

Based on Guilliams et al., Cell 2022. See [Input Data](#input-data) for download instructions.

---

## Scripts

`run_all.py`                            -- runs the full pipeline in order; stops on first failure

  `config.py`                           -- all paths and tunable parameters (N_HVG, N_KMEANS_CLUSTERS, etc.)
  `data_utils.py`                       -- shared read_10x_mtx / write_10x_mtx helpers used by pipeline scripts

  **Pipeline (run in this order)**
  `select_annotated.py`                 -- subsets raw count matrix to cells present in the annotation file
  `select_highly_variable.py`           -- selects top N_HVG highly variable genes (pearson residuals)
  `select_hepatocytes.py`               -- filters to hepatocyte cells using the annot column
  `cluster_hepatocytes.py`              -- normalizes, computes UMAP, k-means clusters, writes per-cluster matrices and umap CSV
  `hepatocytes_scenic.py`               -- runs pySCENIC (GRN -> cisTarget -> AUCell) on each cluster

  **Plots**
  `liver_atlas_umap.py`                 -- UMAP of all annotated cells colored by cell type
  `highly_variable_umap.py`             -- UMAP of HVG-filtered hepatocytes colored by cell type
  `regulon_network.py`                  -- static matplotlib network of top TFs and targets per cluster
  `regulon_heatmap.py`                  -- seaborn clustermap of regulon AUC activity across clusters
  `make_regulon_target_cell_heatmap.py` -- per-cell heatmap of regulon target gene expression, grouped by cluster
  `visualize_networks.py`               -- interactive pyvis HTML network per cluster

---

## Data Directory

`data/source/`
  `annot_humanCD45neg.csv`                -- original cell annotation file; columns: cell barcode, UMAP coords, cell type (annot), patient, diet, sample
  `countTable_human/`                     -- raw 10x count matrix for all human CD45-negative cells (32,738 genes)
  `scenic_databases/`                     -- SCENIC reference files for human (hg38)
    `allTFs_hg38.txt`                     -- list of known human transcription factors
    `hg38_10kbp_up_10kbp_down_...feather` -- cisTarget gene ranking database (10kb window)
    `hg38_500bp_up_100bp_down_...feather` -- cisTarget gene ranking database (500bp window)
    `motifs-v10nr_clust-...tbl`           -- motif annotation table mapping motifs to TFs

`data/intermediate/`
  `selected_cells/`                       -- count matrix subsetted to annotated cells only (select_annotated.py)
  `hv_counts_cells/`                      -- count matrix filtered to top 4000 highly variable genes (select_highly_variable.py)
  `hepatocytes/`                          -- count matrix filtered to hepatocyte cells only (select_hepatocytes.py)
  `hepatocytes_clusters/cluster_{0,1,2}/` -- per-cluster count matrices from k-means clustering (cluster_hepatocytes.py)
  `hepatocytes_scenic/cluster_{0,1,2}/`   -- SCENIC output per cluster: adjacencies.csv, regulons.csv, motifs.csv, auc_matrix.csv

`data/output/`
  `liver_atlas_umap.png`                     -- UMAP of all annotated cells colored by cell type
  `highly_variable_umap.png`                 -- UMAP of hepatocytes after HVG selection, colored by cell type
  `hepatocyte_umap.png`                      -- UMAP of hepatocytes colored by k-means cluster
  `umap_hepatocytes_kmeans.csv`              -- UMAP coordinates + cluster assignments + cell barcodes for hepatocytes
  `network_cluster_{0,1,2}.html`             -- interactive pyvis regulon network per cluster
  `regulon_network_cluster_{0,1,2}.png`      -- static matplotlib regulon network per cluster (top TFs by target count)
  `regulon_target_cell_heatmap.png/.pdf`     -- heatmap of regulon target gene expression per cell, grouped by cluster
  `regulon_target_gene_block_assignment.csv` -- per-gene cluster block assignment from the heatmap

---

## Input Data

**Study data** (`annot_humanCD45neg.csv`, `countTable_human/`):
  Guilliams et al., Cell 2022
  "Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches"
  Data available via the scripts and supplementary materials at:
    https://github.com/saeyslab/LiverCellAtlas (scripts_GuilliamsEtAll_Cell2022)
  Raw count data deposited on GEO -- see the paper's Data Availability section for the accession number.

**SCENIC databases** (`data/source/scenic_databases/`):
  Download from the Aerts lab cisTarget resource:
    https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/
  Files needed:
    `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
    `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
    `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
    `allTFs_hg38.txt`

Your source directory should look like:
```
data/source/
data/source/annot_humanCD45neg.csv
data/source/countTable_human/
data/source/countTable_human/barcodes.tsv.gz
data/source/countTable_human/features.tsv.gz
data/source/countTable_human/matrix.mtx.gz
data/source/scenic_databases/
data/source/scenic_databases/allTFs_hg38.txt
data/source/scenic_databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
data/source/scenic_databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
data/source/scenic_databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

---

## Python Environment

Python 3.9.25

| Package      | Version   |
|--------------|-----------|
| anndata      | 0.10.9    |
| scanpy       | 1.10.3    |
| pyscenic     | 0.12.1    |
| ctxcore      | 0.2.0     |
| scipy        | 1.13.1    |
| numpy        | 2.0.2     |
| pandas       | 2.3.3     |
| scikit-learn | 1.6.1     |
| umap-learn   | 0.5.12    |
| matplotlib   | 3.9.4     |
| seaborn      | 0.13.2    |
| networkx     | 3.2.1     |
| pyvis        | 0.3.2     |
| dask         | 2024.8.0  |
| numba        | 0.60.0    |
| arboreto     | 0.1.6     |

```
pip install anndata==0.10.9 scanpy==1.10.3 pyscenic==0.12.1 ctxcore==0.2.0 scipy==1.13.1 numpy==2.0.2 pandas==2.3.3 scikit-learn==1.6.1 umap-learn==0.5.12 matplotlib==3.9.4 seaborn==0.13.2 networkx==3.2.1 pyvis==0.3.2 "dask[dataframe]==2024.8.0" numba==0.60.0 arboreto==0.1.6
```

Full environment: `project-env/` (run `pip freeze` for complete list)
