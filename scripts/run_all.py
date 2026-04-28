import subprocess
import sys
import os

PYTHON = sys.executable
HERE   = os.path.dirname(os.path.abspath(__file__))

def run(script):
    print(f"\n{'='*60}")
    print(f"  {script}")
    print(f"{'='*60}")
    result = subprocess.run([PYTHON, os.path.join(HERE, script)])
    if result.returncode != 0:
        print(f"\nERROR: {script} failed (exit {result.returncode}). Stopping.")
        sys.exit(result.returncode)

# ── 1. Subset raw counts to annotated cells
run("select_annotated.py")

# ── 2. Select top highly variable genes
run("select_highly_variable.py")

# ── 3. Filter to hepatocytes only
run("select_hepatocytes.py")

# ── 4. Cluster hepatocytes, write per-cluster matrices and UMAP CSV
run("cluster_hepatocytes.py")

# ── 5. Rank marker genes per cluster (Wilcoxon)
run("rank_genes.py")

# ── 6. Run SCENIC regulon inference per cluster
run("hepatocytes_scenic.py")

# ── 6. Plots
run("liver_atlas_umap.py")
run("highly_variable_umap.py")
run("regulon_network.py")
run("regulon_heatmap.py")
run("make_regulon_target_cell_heatmap.py")
run("visualize_networks.py")

print("\nAll steps complete.")
