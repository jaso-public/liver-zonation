import os
import scipy.io
import scipy.sparse
import pandas as pd
import anndata as ad


def _find(directory, base):
    for name in [f'{base}.gz', base]:
        path = os.path.join(directory, name)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(f"{base}(.gz) not found in {directory}")


def read_10x_mtx(directory):
    matrix   = scipy.io.mmread(_find(directory, 'matrix.mtx')).T.tocsr()
    barcodes = pd.read_csv(_find(directory, 'barcodes.tsv'), header=None, sep='\t')[0]
    features = pd.read_csv(_find(directory, 'features.tsv'), header=None, sep='\t')[0]

    adata = ad.AnnData(X=matrix)
    adata.obs_names = barcodes.values
    adata.var_names = features.values
    print(f"Read:    {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def write_10x_mtx(adata, directory):
    os.makedirs(directory, exist_ok=True)
    pd.Series(adata.obs_names).to_csv(os.path.join(directory, 'barcodes.tsv'), sep='\t', header=False, index=False)
    pd.Series(adata.var_names).to_csv(os.path.join(directory, 'features.tsv'), sep='\t', header=False, index=False)
    print(f"Writing matrix ({adata.n_obs} cells x {adata.n_vars} genes)...")
    scipy.io.mmwrite(os.path.join(directory, 'matrix.mtx'), scipy.sparse.csr_matrix(adata.X).T)
    print(f"Written: {adata.n_obs} cells x {adata.n_vars} genes -> {directory}")
