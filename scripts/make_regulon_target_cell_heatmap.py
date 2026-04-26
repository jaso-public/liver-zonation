"""
make_regulon_target_cell_heatmap.py
-----------------------------------
Cell-level heatmap of regulon target-gene expression.

Rows    = union of regulon target genes across the three cluster-specific
          SCENIC runs, grouped into blocks by the cluster where the gene
          is most highly expressed (among clusters whose regulons list it
          as a target).
Columns = individual hepatocyte cells, ordered by cluster
          (cluster_0 | cluster_1 | cluster_2) with thin dividers between
          cluster groups.
Values  = log1p(CP10K) expression per cell.
          Optional row z-score (default ON) so relative differences pop.

Inputs
  - 10x count matrix:   matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz
  - SCENIC regulons:    <scenic-dir>/cluster_*/regulons.csv
  - UMAP cluster CSV:   UMAP_1, UMAP_2, kmeans_cluster, cell

Outputs (under --out-dir)
  - regulon_target_cell_heatmap.png
  - regulon_target_cell_heatmap.pdf
  - regulon_target_gene_block_assignment.csv (written alongside)

Usage
  python make_regulon_target_cell_heatmap.py
  python make_regulon_target_cell_heatmap.py --subsample-per-cluster 1000
  python make_regulon_target_cell_heatmap.py --raw-values
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import config as _config


# ---------------------------------------------------------------- utils


def _find(directory: Path, base: str) -> Path:
    for name in [f"{base}.gz", base]:
        p = directory / name
        if p.exists():
            return p
    raise FileNotFoundError(f"{base}(.gz) not found in {directory}")





def stream_mtx(mtx_path: Path, chunksize: int = 10_000_000):
    reader = pd.read_csv(
        mtx_path,
        sep=" ",
        header=None,
        names=["row", "col", "val"],
        dtype={"row": np.int32, "col": np.int32, "val": np.int32},
        comment="%",
        skiprows=3,
        chunksize=chunksize,
        engine="c",
        compression="infer",
    )
    for chunk in reader:
        yield (
            chunk["row"].to_numpy() - 1,
            chunk["col"].to_numpy() - 1,
            chunk["val"].to_numpy(),
        )


# ---------------------------------------------------------------- main logic

def load_regulon_targets(scenic_dir: Path, clusters: list[str]) -> dict[str, set[str]]:
    gene_to_clusters: dict[str, set[str]] = {}
    for c in clusters:
        df = pd.read_csv(scenic_dir / c / "regulons.csv")
        for g in df["target"].unique():
            gene_to_clusters.setdefault(g, set()).add(c)
    return gene_to_clusters


def build_cell_expression(
    mtx_path: Path,
    n_cells: int,
    target_gene_rows: np.ndarray,       # original gene indices of targets
    hepatocyte_local_idx: np.ndarray,   # len n_cells; -1 if not in plot, else 0..K-1
    chunksize: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Single-pass build of
        expr    : (n_targets, n_hepatocytes) float32, log1p(CP10K)
        (also returns cell_totals for diagnostics)
    """
    n_hepato = int(hepatocyte_local_idx.max()) + 1
    n_targets = len(target_gene_rows)

    target_lookup = np.full(n_cells := len(hepatocyte_local_idx), -1,
                            dtype=np.int32)  # unused; just a trick to keep lints quiet
    # Gene-row lookup: original_gene_idx -> target-row position (or -1)
    max_gene_row = int(target_gene_rows.max()) + 1
    target_row_of = np.full(max_gene_row + 1, -1, dtype=np.int32)
    # But target_gene_rows can come from a larger universe; use full gene count.
    # (We pad to cover the biggest index.)
    # Note: the caller guarantees target_gene_rows are valid indices into features.

    # We'll accumulate counts per (target_row, hepatocyte_local_col) as float32.
    # Memory: 2095 * 15481 * 4 bytes ~= 130 MB. Manageable.
    expr = np.zeros((n_targets, n_hepato), dtype=np.float32)
    cell_totals = np.zeros(n_hepato, dtype=np.float64)

    # Build the target-row lookup correctly (over the full gene-index range used by the mtx)
    # We need a lookup indexed by the original gene index; use the maximum observed index.
    # The safest is to size it to cover the largest index that could appear in the file.
    # The caller sets us up with target_gene_rows in [0, n_genes), so max+1 is enough.
    target_row_of = np.full(int(target_gene_rows.max()) + 2, -1, dtype=np.int32)
    target_row_of[target_gene_rows] = np.arange(n_targets, dtype=np.int32)

    print("streaming mtx: accumulating raw counts for target rows & per-cell totals")
    for row, col, val in stream_mtx(mtx_path, chunksize=chunksize):
        local = hepatocyte_local_idx[col]
        cell_mask = local >= 0
        if not cell_mask.any():
            continue
        local_c = local[cell_mask]
        val_c = val[cell_mask].astype(np.float64)
        # total UMI per hepatocyte (all genes)
        cell_totals += np.bincount(local_c, weights=val_c, minlength=n_hepato)
        # raw counts for target genes only
        r = row[cell_mask]
        # Guard against indices beyond the lookup (shouldn't happen, but safe).
        in_range = r < target_row_of.shape[0]
        if not in_range.all():
            r = r[in_range]
            local_c = local_c[in_range]
            val_c = val_c[in_range]
        tgt = target_row_of[r]
        tmask = tgt >= 0
        if not tmask.any():
            continue
        # scatter into expr: can't use bincount easily for 2-D;
        # use np.add.at, which is fine given only ~0.3% of entries remain.
        np.add.at(expr, (tgt[tmask], local_c[tmask]), val_c[tmask].astype(np.float32))
    print("stream done")

    # CP10K normalization then log1p, done column-wise in-place to save memory.
    scale = np.zeros_like(cell_totals)
    nz = cell_totals > 0
    scale[nz] = 1e4 / cell_totals[nz]
    expr *= scale.astype(np.float32)[np.newaxis, :]
    np.log1p(expr, out=expr)

    return expr, cell_totals


# ---------------------------------------------------------------- plotting

def plot_cell_heatmap(
    expr: np.ndarray,                  # (n_targets, n_hepato)
    genes: list[str],
    gene_block: pd.Series,             # index=genes, values in clusters
    cell_cluster: np.ndarray,          # len n_hepato, cluster id 0..K-1
    clusters: list[str],
    out_dir: Path,
    zscore: bool,
) -> None:
    # Sort rows by block then by each block's mean expression within that block
    per_cluster_mean = np.vstack(
        [expr[:, cell_cluster == i].mean(axis=1) for i in range(len(clusters))]
    ).T                                # (n_targets, n_clusters)

    gene_order: list[int] = []
    for i, c in enumerate(clusters):
        in_block = np.where(gene_block.to_numpy() == c)[0]
        order = in_block[np.argsort(-per_cluster_mean[in_block, i])]
        gene_order.extend(order.tolist())
    gene_order = np.array(gene_order, dtype=np.int64)
    expr = expr[gene_order]
    genes = [genes[i] for i in gene_order]
    gene_block = gene_block.iloc[gene_order]

    # Sort columns by cluster; within cluster, order by PC1 so similar cells sit next to each other
    col_order_list: list[np.ndarray] = []
    divider_positions: list[int] = []
    running = 0
    for i in range(len(clusters)):
        idx = np.where(cell_cluster == i)[0]
        # cheap intra-cluster ordering: sort by total expression of target genes
        within_score = expr[:, idx].sum(axis=0)
        idx = idx[np.argsort(within_score)]
        col_order_list.append(idx)
        running += len(idx)
        divider_positions.append(running)
    col_order = np.concatenate(col_order_list)
    expr = expr[:, col_order]
    cell_cluster = cell_cluster[col_order]

    # Prepare display matrix
    if zscore:
        mu = expr.mean(axis=1, keepdims=True)
        sd = expr.std(axis=1, keepdims=True)
        sd[sd == 0] = np.nan
        display = (expr - mu) / sd
        vmin, vmax = -2.5, 2.5
        cbar_label = "Row z-score of log1p(CP10K)"
    else:
        display = expr
        vmin = 0.0
        vmax = float(np.nanpercentile(display, 99))
        cbar_label = "log1p(CP10K) expression"

    cmap = LinearSegmentedColormap.from_list(
        "scenic", ["#08306b", "#f7fbff", "#b30000"], N=256
    )
    cmap.set_bad("#dddddd")

    # Block colors
    _palette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                "#a65628", "#f781bf", "#999999", "#ffff33", "#a6cee3"]
    block_colors = {c: _palette[i % len(_palette)] for i, c in enumerate(clusters)}
    gene_block_rgb = np.array([plt.matplotlib.colors.to_rgb(block_colors[b])
                               for b in gene_block])
    cell_cluster_rgb = np.array([plt.matplotlib.colors.to_rgb(
        block_colors[clusters[i]]) for i in cell_cluster])

    # Figure with 4 axes: top cluster bar, left gene-block bar, main heatmap, colorbar
    n_rows = expr.shape[0]
    n_cols = expr.shape[1]
    fig_h = max(10, min(0.02 * n_rows + 4, 40))
    fig = plt.figure(figsize=(14, fig_h))
    gs = fig.add_gridspec(
        2, 3,
        width_ratios=[0.015, 1, 0.02],
        height_ratios=[0.012, 1],
        wspace=0.01, hspace=0.01,
    )
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_main = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 2])

    # top bar (cells colored by cluster)
    ax_top.imshow(cell_cluster_rgb.reshape(1, -1, 3), aspect="auto")
    ax_top.set_xticks([]); ax_top.set_yticks([])
    for i, c in enumerate(clusters):
        idx = np.where(cell_cluster == i)[0]
        if len(idx) == 0:
            continue
        mid = int(idx.mean())
        ax_top.text(mid, -0.6, c.replace("_", " ").title(),
                    ha="center", va="bottom", fontsize=10, color=block_colors[c])

    # left bar (genes colored by assigned block)
    ax_left.imshow(gene_block_rgb.reshape(-1, 1, 3), aspect="auto")
    ax_left.set_xticks([]); ax_left.set_yticks([])
    ax_left.set_ylabel("Gene block", fontsize=9)

    # main heatmap
    im = ax_main.imshow(
        display, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest",
    )
    ax_main.set_xticks([]); ax_main.set_yticks([])
    ax_main.set_xlabel(f"{n_cols} cells (grouped by cluster)")
    ax_main.set_ylabel(f"{n_rows} regulon target genes")
    ax_main.set_title("Regulon target-gene expression per cell\n"
                      "cells grouped by cluster",
                      fontsize=12, pad=12)

    # white dashed dividers between cell-clusters
    for x in divider_positions[:-1]:
        ax_main.axvline(x - 0.5, color="white", linewidth=1.2)

    # colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label(cbar_label, fontsize=10)

    # legend
    handles = [plt.Line2D([0], [0], marker="s", linestyle="",
                          markerfacecolor=block_colors[c], markeredgecolor="k",
                          markersize=9, label=c.replace("_", " ").title())
               for c in clusters]
    ax_main.legend(handles=handles, title="Cluster",
                   bbox_to_anchor=(1.10, 1.0), loc="upper left",
                   frameon=False, fontsize=9, title_fontsize=9)

    png = out_dir / "regulon_target_cell_heatmap.png"
    pdf = out_dir / "regulon_target_cell_heatmap.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {pdf}")


# ---------------------------------------------------------------- driver

count_dir  = Path(_config.hepatocytes_count_dir)
scenic_dir = Path(_config.hepatocytes_scenic_dir)
umap_csv   = Path(_config.OUTPUT_DIR) / "umap_hepatocytes_kmeans.csv"
out_dir    = Path(_config.OUTPUT_DIR)
clusters   = [f"cluster_{i}" for i in range(_config.N_KMEANS_CLUSTERS)]

out_dir.mkdir(parents=True, exist_ok=True)
# ---- metadata ---------------------------------------------------------
print("Loading barcodes / features / UMAP clusters ...")
barcodes = pd.read_csv(_find(count_dir, "barcodes.tsv"), header=None, sep='\t')[0].tolist()
features = pd.read_csv(_find(count_dir, "features.tsv"), header=None, sep='\t')[0].tolist()
umap = pd.read_csv(umap_csv)
print(f"  {len(barcodes)} cells, {len(features)} genes, "
      f"{len(umap)} hepatocytes in umap")

barcode_idx = {b: i for i, b in enumerate(barcodes)}
gene_idx = {g: i for i, g in enumerate(features)}

umap = umap[umap["cell"].isin(barcode_idx)].copy()
umap = umap.sort_values(["kmeans_cluster", "cell"]).reset_index(drop=True)

# Build hepatocyte_local_idx: len = total cells in mtx
hepatocyte_local_idx = np.full(len(barcodes), -1, dtype=np.int32)
cell_cluster_local: list[int] = []
for local_i, (cell, clust) in enumerate(zip(umap["cell"], umap["kmeans_cluster"])):
    hepatocyte_local_idx[barcode_idx[cell]] = local_i
    cell_cluster_local.append(int(clust))
cell_cluster = np.array(cell_cluster_local, dtype=np.int32)
print(f"  plotting {len(umap)} cells")
for c in clusters:
    k = int(c.split("_")[1])
    print(f"    {c}: {(cell_cluster == k).sum()} cells")

# ---- target genes -----------------------------------------------------
print("Loading regulon targets ...")
gene_to_clusters = load_regulon_targets(scenic_dir, clusters)
target_genes = sorted(g for g in gene_to_clusters if g in gene_idx)
target_rows = np.array([gene_idx[g] for g in target_genes], dtype=np.int32)
print(f"  {len(target_genes)} target genes found in the count matrix")

# ---- per-cell expression ---------------------------------------------
print("Streaming count matrix ...")
expr, _ = build_cell_expression(
    mtx_path=_find(count_dir, "matrix.mtx"),
    n_cells=len(barcodes),
    target_gene_rows=target_rows,
    hepatocyte_local_idx=hepatocyte_local_idx,
    chunksize=20_000_000,
)

# ---- per-gene block assignment ----------------------------------------
per_cluster_mean = np.vstack(
    [expr[:, cell_cluster == int(c.split("_")[1])].mean(axis=1)
     for c in clusters]
).T
name_to_idx = {c: i for i, c in enumerate(clusters)}
block_assignment = []
for g, row in zip(target_genes, per_cluster_mean):
    allowed_idx = [name_to_idx[c] for c in gene_to_clusters[g] if c in name_to_idx]
    masked = np.full_like(row, -np.inf, dtype=float)
    masked[allowed_idx] = row[allowed_idx]
    block_assignment.append(clusters[int(np.argmax(masked))])
block = pd.Series(block_assignment, index=target_genes, name="block")

# Filter 1: keep only genes unique to one cluster
unique_mask = np.array([len(gene_to_clusters[g]) == 1 for g in target_genes])
# Filter 2: keep top N_HEATMAP_GENES by variance across cells
gene_var = expr.var(axis=1)
top_var_mask = np.zeros(len(target_genes), dtype=bool)
top_var_mask[np.argsort(gene_var)[-_config.N_HEATMAP_GENES:]] = True

keep = unique_mask & top_var_mask
expr         = expr[keep]
target_genes = [g for g, k in zip(target_genes, keep) if k]
block        = block[keep]
print(f"  {keep.sum()} genes after filtering (cluster-unique + top {_config.N_HEATMAP_GENES} by variance)")

block.to_frame().assign(
    membership=[",".join(sorted(gene_to_clusters[g])) for g in block.index]
).to_csv(out_dir / "regulon_target_gene_block_assignment.csv")

# ---- plot -------------------------------------------------------------
plot_cell_heatmap(
    expr=expr,
    genes=target_genes,
    gene_block=block,
    cell_cluster=cell_cluster,
    clusters=clusters,
    out_dir=out_dir,
    zscore=True,
)
print("all done")
