import config

import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

os.makedirs(config.OUTPUT_DIR, exist_ok=True)

for cluster_id in range(config.N_KMEANS_CLUSTERS):
    regulon_path = os.path.join(config.hepatocytes_scenic_dir, f"cluster_{cluster_id}", "regulons.csv")
    if not os.path.exists(regulon_path):
        print(f"Cluster {cluster_id}: regulons.csv not found, skipping")
        continue

    regulons = pd.read_csv(regulon_path)
    regulons['tf_name'] = regulons['regulon'].str.split('(').str[0].str.strip()

    # Pick top N_TOP_TFS by number of targets
    top_tfs = (regulons.groupby('tf_name')['target']
               .count()
               .nlargest(config.N_TOP_TFS)
               .index.tolist())
    edges = (regulons[regulons['tf_name'].isin(top_tfs)]
             .groupby('tf_name', group_keys=False)
             .head(config.N_TOP_TARGETS))

    # Build graph
    G = nx.DiGraph()
    for _, row in edges.iterrows():
        G.add_edge(row['tf_name'], row['target'])

    tf_nodes     = [n for n in G.nodes if n in top_tfs]
    target_nodes = [n for n in G.nodes if n not in top_tfs]
    print(f"Cluster {cluster_id}: {len(tf_nodes)} TFs, {len(target_nodes)} targets, {G.number_of_edges()} edges")

    pos = nx.spring_layout(G, seed=42, k=0.5)

    plt.figure(figsize=(14, 10))
    plt.rcParams['savefig.bbox'] = 'tight'
    nx.draw_networkx_nodes(G, pos, nodelist=tf_nodes,     node_color='steelblue', node_size=600)
    nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_color='lightgrey',  node_size=80)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in tf_nodes},     font_size=8, font_weight='bold')
    nx.draw_networkx_labels(G, pos, labels={n: n for n in target_nodes}, font_size=6, font_color='#444444')
    nx.draw_networkx_edges(G, pos, alpha=0.2, arrows=True, arrowsize=8, edge_color='grey')

    plt.legend(handles=[
        mpatches.Patch(color='steelblue', label='TF'),
        mpatches.Patch(color='lightgrey', label='Target gene'),
    ], loc='upper left')
    plt.title(f"Regulon network — cluster {cluster_id} (top {config.N_TOP_TFS} TFs by target count)", fontsize=13)
    plt.axis('off')

    out_path = os.path.join(config.OUTPUT_DIR, f"regulon_network_cluster_{cluster_id}")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Saved -> {out_path}")
