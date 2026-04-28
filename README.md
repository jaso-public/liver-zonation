# Hepatocyte Regulon Analysis

Liver atlas single-cell pipeline: HVG selection → hepatocyte clustering → pySCENIC regulon inference → visualization.

Based on Guilliams et al., Cell 2022. See [Input Data](#input-data) for download instructions.

---

## Getting Started

**Step 1 — Clone the repository**
```
git clone https://github.com/jaso-public/liver-zonation.git
cd liver-zonation
```

**Step 2 — Create a Python virtual environment**
```
python3.9 -m venv lz-env
```

**Step 3 — Activate the environment**
```
source lz-env/bin/activate
```

**Step 4 — Install dependencies**
```
pip install anndata==0.10.9 scanpy==1.10.3 pyscenic==0.12.1 ctxcore==0.2.0 scipy==1.13.1 numpy==2.0.2 pandas==2.3.3 scikit-learn==1.6.1 umap-learn==0.5.12 matplotlib==3.9.4 seaborn==0.13.2 networkx==3.2.1 pyvis==0.3.2 "dask[dataframe]==2024.8.0" numba==0.60.0 arboreto==0.1.6
```

**Step 5 — Download the source data**

See [Input Data](#input-data) for instructions on downloading the study count matrix and SCENIC databases into `data/source/`.

**Step 6 — Run the pipeline**
```
python scripts/run_all.py
```

---

## Scripts

- `run_all.py` — runs the full pipeline in order; stops on first failure
- `config.py` — all paths and tunable parameters (`N_HVG`, `N_KMEANS_CLUSTERS`, etc.)
- `data_utils.py` — shared `read_10x_mtx` / `write_10x_mtx` helpers used across the pipeline

### Pipeline (run in this order)

| Script | Description |
|---|---|
| `select_annotated.py` | Subsets raw count matrix to cells present in the annotation file |
| `select_highly_variable.py` | Selects top `N_HVG` highly variable genes (Pearson residuals) |
| `select_hepatocytes.py` | Filters to hepatocyte cells using the `annot` column |
| `cluster_hepatocytes.py` | Normalizes, computes UMAP, k-means clusters, writes per-cluster matrices and UMAP CSV |
| `rank_genes.py` | Wilcoxon rank-sum test per cluster; outputs ranked marker gene CSV |
| `hepatocytes_scenic.py` | Runs pySCENIC (GRN → cisTarget → AUCell) on each cluster |

### Plots

| Script | Output |
|---|---|
| `plot_marker_genes.py` | Dot plot and log2FC bar charts for top marker genes per cluster |
| `liver_atlas_umap.py` | UMAP of all annotated cells colored by cell type |
| `highly_variable_umap.py` | UMAP of HVG-filtered hepatocytes colored by cell type |
| `regulon_network.py` | Static matplotlib network of top TFs and targets per cluster |
| `regulon_heatmap.py` | Seaborn clustermap of regulon AUC activity across clusters |
| `make_regulon_target_cell_heatmap.py` | Per-cell heatmap of regulon target gene expression, grouped by cluster |
| `visualize_networks.py` | Interactive pyvis HTML network per cluster |

---

## Data Directory

### `data/source/`

| Path | Description |
|---|---|
| `annot_humanCD45neg.csv` | Original cell annotation file; columns: cell barcode, UMAP coords, cell type, patient, diet, sample |
| `countTable_human/` | Raw 10x count matrix for all human CD45-negative cells (32,738 genes) |
| `scenic_databases/allTFs_hg38.txt` | List of known human transcription factors |
| `scenic_databases/hg38_10kbp_...feather` | cisTarget gene ranking database (10kb window) |
| `scenic_databases/hg38_500bp_...feather` | cisTarget gene ranking database (500bp window) |
| `scenic_databases/motifs-v10nr_...tbl` | Motif annotation table mapping motifs to TFs |

### `data/intermediate/`

| Path | Description |
|---|---|
| `selected_cells/` | Count matrix subsetted to annotated cells only |
| `hv_counts_cells/` | Count matrix filtered to top 4000 highly variable genes |
| `hepatocytes/` | Count matrix filtered to hepatocyte cells only |
| `hepatocytes_clusters/cluster_{0,1,2}/` | Per-cluster count matrices from k-means clustering |
| `hepatocytes_scenic/cluster_{0,1,2}/` | SCENIC output: `adjacencies.csv`, `regulons.csv`, `motifs.csv`, `auc_matrix.csv` |

### `data/output/`

| File | Description |
|---|---|
| `liver_atlas_umap.png` | UMAP of all annotated cells colored by cell type |
| `highly_variable_umap.png` | UMAP of hepatocytes after HVG selection, colored by cell type |
| `hepatocyte_umap.png` | UMAP of hepatocytes colored by k-means cluster |
| `umap_hepatocytes_kmeans.csv` | UMAP coordinates, cluster assignments, and cell barcodes for hepatocytes |
| `network_cluster_{0,1,2}.html` | Interactive pyvis regulon network per cluster |
| `regulon_network_cluster_{0,1,2}.png` | Static regulon network per cluster (top TFs by target count) |
| `regulon_target_cell_heatmap.png/.pdf` | Per-cell heatmap of regulon target gene expression, grouped by cluster |
| `regulon_target_gene_block_assignment.csv` | Per-gene cluster block assignment from the heatmap |

---

## Input Data

### Study data (`annot_humanCD45neg.csv`, `countTable_human/`)

Guilliams et al., Cell 2022 — *"Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches"*

- Scripts and supplementary data: (https://www.livercellatlas.org/download.php)
- Raw count data deposited on GEO — see the paper's Data Availability section for the accession number.

### SCENIC databases (`data/source/scenic_databases/`)

Download from the Aerts lab cisTarget resource:
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/

Files needed:
- `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
- `allTFs_hg38.txt`

Your `data/source/` directory should look like:

```
data/source/
├── annot_humanCD45neg.csv
├── countTable_human/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── scenic_databases/
    ├── allTFs_hg38.txt
    ├── hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
    ├── hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
    └── motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

---

## Python Environment

Python 3.9.25

| Package | Version |
|---|---|
| anndata | 0.10.9 |
| scanpy | 1.10.3 |
| pyscenic | 0.12.1 |
| ctxcore | 0.2.0 |
| scipy | 1.13.1 |
| numpy | 2.0.2 |
| pandas | 2.3.3 |
| scikit-learn | 1.6.1 |
| umap-learn | 0.5.12 |
| matplotlib | 3.9.4 |
| seaborn | 0.13.2 |
| networkx | 3.2.1 |
| pyvis | 0.3.2 |
| dask | 2024.8.0 |
| numba | 0.60.0 |
| arboreto | 0.1.6 |

```
pip install anndata==0.10.9 scanpy==1.10.3 pyscenic==0.12.1 ctxcore==0.2.0 scipy==1.13.1 numpy==2.0.2 pandas==2.3.3 scikit-learn==1.6.1 umap-learn==0.5.12 matplotlib==3.9.4 seaborn==0.13.2 networkx==3.2.1 pyvis==0.3.2 "dask[dataframe]==2024.8.0" numba==0.60.0 arboreto==0.1.6
```

Full environment: `project-env/` (run `pip freeze` for complete list)
