import pandas as pd
import networkx as nx

data = {'1': ["M830T1198 M374T1008 M552T1006 Solyc02g037500.1 Solyc08g080620.1 Solyc10g075150.1 Solyc01g088810.2"],
        '2': ["M243T818 M283T818 Solyc02g037500.1 Solyc02g067520.2 Solyc08g080620.1 Solyc06g060110.2 Solyc01g087570.2"],
        '3': ["M830T1198 M374T1008 Solyc02g037500.1 Solyc01g110110.2 Solyc08g080620.1"],
        '4': ["M552T1006 M830T1198 M328T929 Solyc02g037500.1 M352T921 Solyc08g080620.1 Solyc10g051230.1 Solyc03g113580.1 Solyc01g059880.2"],
        '5': ["M830T1198 M374T1008 Solyc02g037500.1 Solyc08g007110.2 Solyc08g080620.1"],
        '6': ["M302T590 M376T548 Solyc10g006870.1 Solyc04g011750.2 Solyc08g076310.2 Solyc01g095610.2 Solyc01g079580.2"]}

data_list = [{'ID': key, 'Members': value[0]} for key, value in data.items()]
df = pd.DataFrame(data_list)
df.set_index('ID')

# Split the 'components' column into a list of genes and metabolites
#df['Members_list'] = df['Members'].apply(lambda x: x.split())
new_df = df.assign(Members_list=df['Members'].apply(lambda x: x.split()))

# Create a graph to represent connections between genes and metabolites
G = nx.Graph()

# Add nodes to the graph
G.add_nodes_from(new_df['ID'])

    # Add edges between nodes in different clusters if they share at least one element
    #for i, (_, row) in enumerate(new_df.iterrows()):
for i, (_, row) in enumerate(new_df.iterrows()):
    cluster_nodes = row['Members_list']
    for _, other_row in new_df.iloc[i + 1:].iterrows():
        other_cluster_nodes = other_row['Members_list']
        if any(node in other_cluster_nodes for node in cluster_nodes):
            G.add_edge(row['ID'], other_row['ID'])

# Find connected components (clusters) in the graph
connected_components = list(nx.connected_components(G))

# Create a mapping from gene/metabolite to cluster ID
node_cluster_mapping = {node: cluster_id for cluster_id, nodes in enumerate(connected_components, start=1) for node in nodes}


# Display the merged clusters with original cluster IDs
merged_clusters = {}
for i, cluster in enumerate(connected_components, start=1):
    original_cluster_ids = {int(node_cluster_mapping[node]) for node in cluster}
    merged_cluster_id = f'MC{i}'
    merged_clusters[merged_cluster_id] = cluster

merged_list = [{'Cluster': key, 'Members': ', '.join(value)} for key, value in merged_clusters.items()]
merged_df = pd.DataFrame(merged_list)
print(merged_df)

# Display the merged clusters with original cluster IDs and row IDs
#for _, merged_cluster in merged_clusters.items():
#    print(f"Merged Cluster {merged_cluster}")
