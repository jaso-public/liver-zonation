"""
Hepatocyte regulon analysis — pySCENIC pipeline
================================================
Prerequisites
-------------
1. Install pySCENIC:
       pip install pyscenic

2. Download SCENIC databases for human (hg38) from:
       https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/

   You need:
     - hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
     - hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
     - motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl   (motif annotation)
     - allTFs_hg38.txt                               (human TF list)

   All four are also bundled in the pySCENIC tutorial data:
       https://github.com/aertslab/pySCENIC#tutorial

Pipeline
--------
  Step 1 — GRNBoost2   : infer TF co-expression modules
  Step 2 — cisTarget   : prune modules to motif-supported regulons
  Step 3 — AUCell      : score each cell for regulon activity
  Step 4 — Export      : per-cluster regulon activity summary
"""

import os
import logging
import warnings
warnings.simplefilter('ignore')
logging.getLogger('pyscenic.transform').setLevel(logging.ERROR)
logging.getLogger('pyscenic.utils').setLevel(logging.ERROR)
logging.getLogger('pyscenic.prune').setLevel(logging.ERROR)

import numpy as np

# pyscenic uses np.object/np.bool/np.int/np.float which were removed in NumPy 1.24
np.object = object
np.bool   = np.bool_
np.int    = np.int_
np.float  = np.float64

import pandas as pd
import scipy.sparse

import config
from data_utils import read_10x_mtx

os.makedirs(config.hepatocytes_scenic_dir, exist_ok=True)

if __name__ == '__main__':
    # dask-expr incompatibility fix: from_delayed receives a generator but calls len() on it
    import dask.dataframe as _ddd
    _orig_from_delayed = _ddd.from_delayed
    _ddd.from_delayed = lambda dfs, **kw: _orig_from_delayed(list(dfs), **kw)

    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.aucell import aucell

    tf_names_all = pd.read_csv(config.tf_list, header=None)[0].tolist()
    dbs = [RankingDatabase(fname=f, name=os.path.basename(f)) for f in config.ranking_dbs]

    for cluster_id in range(config.N_KMEANS_CLUSTERS):
        print(f"\n{'='*60}")
        print(f"  CLUSTER {cluster_id}")
        print(f"{'='*60}")

        cl_dir  = os.path.join(config.hepatocytes_clusters_dir, f"cluster_{cluster_id}")
        out_cl  = os.path.join(config.hepatocytes_scenic_dir, f"cluster_{cluster_id}")
        os.makedirs(out_cl, exist_ok=True)

        # ── Load this cluster's matrix ─────────────────────────────────────────
        adata = read_10x_mtx(cl_dir)
        print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

        # matrices are already normalized+scaled from the pipeline — skip normalize/log1p
        ex_matrix = pd.DataFrame(
            adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X,
            index=adata.obs_names,
            columns=adata.var_names,
        )

        # Select top N_SCENIC_GENES by variance, always keeping TFs
        tf_names = [t for t in tf_names_all if t in ex_matrix.columns]
        top_var_genes = ex_matrix.var(0).nlargest(config.N_SCENIC_GENES).index
        keep_genes = top_var_genes.union(pd.Index(tf_names))
        ex_matrix = ex_matrix[keep_genes]
        tf_names = [t for t in tf_names if t in ex_matrix.columns]
        print(f"  Selected {len(ex_matrix.columns)} genes (top {config.N_SCENIC_GENES} HV + {len(tf_names)} TFs)")


        # ── Step 1: TF-target correlations ────────────────────────────────────
        print("  Step 1: TF-target correlations...")

        # Drop zero-variance genes (would produce NaN correlations)
        gene_std = ex_matrix.std(0)
        ex_matrix = ex_matrix.loc[:, gene_std > 0]
        tf_names  = [t for t in tf_names if t in ex_matrix.columns]
        print(f"    {ex_matrix.shape[1]} genes with variance, {len(tf_names)} TFs retained")

        tf_mat   = ex_matrix[tf_names].values   # float64
        gene_mat = ex_matrix.values

        tf_std   = tf_mat.std(0)
        gene_std = gene_mat.std(0)
        tf_z     = (tf_mat   - tf_mat.mean(0))   / np.where(tf_std   > 0, tf_std,   1)
        gene_z   = (gene_mat - gene_mat.mean(0)) / np.where(gene_std > 0, gene_std, 1)
        corr     = (tf_z.T @ gene_z) / len(ex_matrix)   # shape: TFs × genes

        # Flatten into adjacency table explicitly
        tfs_col     = np.repeat(tf_names, len(ex_matrix.columns))
        targets_col = np.tile(ex_matrix.columns.tolist(), len(tf_names))
        imp_col     = corr.ravel()

        adjacencies = pd.DataFrame({'TF': tfs_col, 'target': targets_col, 'importance': imp_col})
        adjacencies = adjacencies[adjacencies['TF'] != adjacencies['target']].copy()
        adjacencies = adjacencies.sort_values('importance', ascending=False).reset_index(drop=True)
        adjacencies.to_csv(os.path.join(out_cl, "adjacencies.csv"), index=False)
        print(f"    {len(adjacencies)} TF-target edges (top importance: {adjacencies['importance'].max():.4f})")

        # ── Step 2: cisTarget pruning ──────────────────────────────────────────
        print("  Step 2: cisTarget pruning...")
        modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
        print(f"    {len(modules)} modules")

        if len(modules) == 0:
            print("    WARNING: no modules — skipping this cluster")
            continue

        df_motifs = prune2df(
            rnkdbs=dbs,
            modules=modules,
            motif_annotations_fname=config.motif_annot,
            num_workers=1,
        )
        df_motifs.to_csv(os.path.join(out_cl, "motifs.csv"))

        regulons = df2regulons(df_motifs)
        print(f"    {len(regulons)} regulons")

        # Save regulons as a flat CSV: TF, target_gene
        reg_rows = [{'regulon': r.name, 'target': gene} for r in regulons for gene in r.genes]
        pd.DataFrame(reg_rows).to_csv(os.path.join(out_cl, "regulons.csv"), index=False)
        print(f"    Regulons written -> {out_cl}/regulons.csv")

        # ── Step 3: AUCell ─────────────────────────────────────────────────────
        print("  Step 3: AUCell scoring...")
        auc_mtx = aucell(exp_mtx=ex_matrix, signatures=regulons, num_workers=1)
        auc_mtx.to_csv(os.path.join(out_cl, "auc_matrix.csv"))
        print(f"    AUC matrix: {auc_mtx.shape}")

        # Top 5 regulons by mean activity
        top5 = auc_mtx.mean().nlargest(5).index.tolist()
        print(f"    Top 5 regulons: {top5}")

    print("\nDone.")
    print(f"Results in: {config.hepatocytes_scenic_dir}")
