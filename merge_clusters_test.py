import pandas as pd


# Read the space-separated file into a DataFrame
df = pd.read_csv('jeon_bac_N_10.csv', sep=',')

# Extract the last column and convert it into a dictionary
cluster_dict = {}
for index, row in df.iterrows():
    clusters = row['Members'].split(' ')
    cluster_key = f'cluster{index + 1}'
    cluster_dict[cluster_key] = clusters

print(cluster_dict)
# Create a list of clusters
clusters = [set(cluster_dict[key]) for key in cluster_dict]

# store M* elements
m_elements = {}

# populate the dictionary
for i, cluster in enumerate(clusters):
    m_star = [element for element in cluster if element.startswith('M')]
    for m in m_star:
        if m in m_elements:
            m_elements[m].append(i)
        else:
            m_elements[m] = [i]

# Merge clusters based on shared M* elements
merged_clusters = []
visited_clusters = set()

for m, indices in m_elements.items():
    if len(indices) > 1:
        merged_set = set()
        for index in indices:
            visited_clusters.add(index)
            merged_set.update(clusters[index])
        #if merged_set not in merged_clusters:
        #    merged_clusters.append(merged_set)

        # Check for overlaps with existing merged clusters
        new_merged_set = merged_set
        for existing_cluster in merged_clusters:
            if existing_cluster & new_merged_set:
                existing_cluster.update(new_merged_set)
                visited_clusters.update(indices)
                break
        else:
            if new_merged_set not in merged_clusters:
                merged_clusters.append(new_merged_set)

# Add remaining clusters that don't share M* elements
for i, cluster in enumerate(clusters):
    if i not in visited_clusters:
        merged_clusters.append(cluster)

# Print the merged clusters
#for i, merged_cluster in enumerate(merged_clusters, start=1):
#    print(f'Merged Cluster {i}: {merged_cluster}')

# Create a list to store metabolites and genes for each merged cluster
metabolites_list = []
genes_list = []

# Iterate through merged clusters and separate elements into metabolites and genes
for i, cluster in enumerate(merged_clusters):
    metabolites = [element for element in cluster if element.startswith('M')]
    genes = [element for element in cluster if element.startswith('S')]
    
    # Append the metabolites and genes to their respective lists
    metabolites_list.append(' '.join(metabolites))
    genes_list.append(' '.join(genes))

# Create a DataFrame with ID, Metabolites, and Genes columns
df = pd.DataFrame({
    'ID': [f'Merged_cluster_{i+1}' for i in range(len(merged_clusters))],
    'Metabolites': metabolites_list,
    'Genes': genes_list
})

# Export the DataFrame to a CSV file
df.to_csv('merged_clusters_jeon_bac_N_10.csv', index=False)



# Convert merged_clusters list to a Pandas DataFrame with two columns
#df = pd.DataFrame({
#    'ID': [f'Merged_cluster_{i+1}' for i in range(len(merged_clusters))],
#    'Elements': [' '.join(cluster) for cluster in merged_clusters]
#})

# Export the DataFrame to a CSV file
#df.to_csv('merged_clusters_jeon_bac_N_10.csv', index=False)
