import networkx as nx
import matplotlib.pyplot as plt

# Create a cyclic graph
G = nx.DiGraph()
G.add_edges_from([(1, 2), (2, 3), (3, 1)])

# Visualize the cyclic graph
pos = nx.circular_layout(G)
nx.draw(G, pos, with_labels=True, font_weight='bold', node_size=700, node_color='skyblue', arrowsize=20)
plt.title('Cyclic Graph')
plt.show()

# Detect and remove cycles
def remove_cycles(graph):
    cycles = list(nx.simple_cycles(graph))
    if cycles:
        cycle = cycles[0]
        print("Cycle found:", cycle)
        # If a cycle is found, remove one edge from the cycle
        edge_to_remove = (cycle[0], cycle[1])
        graph.remove_edge(*edge_to_remove)
        print("Edge removed:", edge_to_remove)
    else:
        print("No cycles found.")

# Detect and remove cycles until there are none left
while True:
    remove_cycles(G)
    if not list(nx.simple_cycles(G)):
        break

# Visualize the acyclic graph
pos = nx.circular_layout(G)
nx.draw(G, pos, with_labels=True, font_weight='bold', node_size=700, node_color='lightgreen', arrowsize=20)
plt.title('Acyclic Graph')
plt.show()
