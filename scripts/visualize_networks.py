import config

import pandas as pd
from pyvis.network import Network
import os

clusters = [f"cluster_{i}" for i in range(config.N_KMEANS_CLUSTERS)]

for cluster in clusters:
    csv_path = os.path.join(config.hepatocytes_scenic_dir, cluster, 'regulons.csv')
    output_html = os.path.join(config.OUTPUT_DIR, f'network_{cluster}.html')
    
    # Check if the file actually exists before processing
    if not os.path.exists(csv_path):
        print(f"Skipping {cluster}: File not found at {csv_path}")
        continue

    print(f"Processing {cluster}...")
    
    # 1. Load data
    df = pd.read_csv(csv_path)
    
    # 2. Setup Network (optimized for a clear view)
    net = Network(
        height='800px', 
        width='100%', 
        directed=True, 
        notebook=False, # Set to False since we are saving files to specific paths
        cdn_resources='remote',
        bgcolor='#ffffff', 
        font_color='black'
    )

    # 3. Add data (limiting to top 200 for performance, remove .head() for all)
    subset = df.head(400)
    
    for _, row in subset.iterrows():
        tf = row['regulon']
        target = row['target']
        
        # Add TF node (Diamond shape to distinguish from genes)
        net.add_node(tf, label=tf, color='#e74c3c', title="Transcription Factor", shape='diamond')
        # Add Target node
        net.add_node(target, label=target, color='#3498db', title="Target Gene")
        # Add Edge
        net.add_edge(tf, target, color='#bdc3c7')

    # 4. Save the file directly into the cluster directory
    net.write_html(output_html)
    print(f"Successfully saved: {output_html}")

print("\nAll visualizations complete!")

