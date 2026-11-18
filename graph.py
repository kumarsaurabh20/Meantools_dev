import networkx as nx
import matplotlib.pyplot as plt

# Define clusters
cluster1 = {'M2', 'S1', 'M3', 'M4', 'S3', 'S4', 'S2', 'S5', 'S26', 'S7', 'M7', 'M6'}
cluster2 = {'M2', 'S1', 'M3', 'M41', 'S3', 'S41', 'S21', 'S51', 'S26', 'S71', 'M71', 'M61'}
cluster3 = {'M2', 'S1', 'M3', 'M24', 'S24', 'S25', 'S221', 'S251', 'S216'}

# Create graphs for each cluster
G1 = nx.Graph()
G1.add_edges_from([(node1, node2) for node1 in cluster1 for node2 in cluster1 if node1 != node2])

G2 = nx.Graph()
G2.add_edges_from([(node1, node2) for node1 in cluster2 for node2 in cluster2 if node1 != node2])

G3 = nx.Graph()
G3.add_edges_from([(node1, node2) for node1 in cluster3 for node2 in cluster3 if node1 != node2])

# Assign node colors
node_colors_G1 = ['blue' if node.startswith('M') else 'green' for node in G1.nodes]
node_colors_G2 = ['blue' if node.startswith('M') else 'green' for node in G2.nodes]
node_colors_G3 = ['blue' if node.startswith('M') else 'green' for node in G3.nodes]

# Plot the graphs with node colors
plt.figure(figsize=(15, 5))

plt.subplot(131)
nx.draw(G1, with_labels=True, font_weight='bold', node_color=node_colors_G1)
plt.title('Cluster 1')

plt.subplot(132)
nx.draw(G2, with_labels=True, font_weight='bold', node_color=node_colors_G2)
plt.title('Cluster 2')

plt.subplot(133)
nx.draw(G3, with_labels=True, font_weight='bold', node_color=node_colors_G3)
plt.title('Cluster 3')

plt.show()
