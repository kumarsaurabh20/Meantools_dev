import os
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm


pca_df = pd.read_csv('Arabidopsis_clusters_DR_25_pca_labelled.csv')
cluster_members_df = pd.read_csv('Ara_clusters_members.csv')
gene_expression_df = pd.read_csv('/Users/singh018/Downloads/NWO_groot_paired_omics_2024/Transcriptomics_Run_2024/TR_DEGs/roots_7dpi/TR_combined_MH_MG_QUANT.csv', index_col=0)  # Replace with actual file
metabolite_abundance_df = pd.read_csv('/Users/singh018/Downloads/NWO_groot_paired_omics_2024/Transcriptomics_Run_2024/TR_DEGs/roots_7dpi/MT_combined_MH_MG_QUANT_filtered_replicates.csv', index_col=0)  # Replace with actual file

#Extract clusters and create folders
for cluster_id in pca_df['Cluster'].unique():
    cluster_df = pca_df[pca_df['Cluster'] == cluster_id]
    
    #Create a folder for each cluster
    cluster_folder = f'cluster_{cluster_id}'
    os.makedirs(cluster_folder, exist_ok=True)

    #Extract genes and metabolites from the cluster members file
    
    for idx, row in tqdm(cluster_df.iterrows(), total=cluster_df.shape[0], desc=f'Processing Cluster {cluster_id}'):
        #for _, row in cluster_df.iterrows():
        current_id = row['ID']
        matching_row = cluster_members_df[cluster_members_df['ID']== current_id]

        if matching_row.empty:
            print(f"No matching members found for ID {current_id} in Cluster {cluster_id}. Skipping.")
            continue

        cluster_members = matching_row['Members'].values[0]
        members_list = cluster_members.split()

        genes_in_cluster = [member for member in members_list if member.startswith("AT")]
        metabolites_in_cluster = [member for member in members_list if member.startswith("Ar")]

        #Plot gene expression pattern for the cluster
        if genes_in_cluster:
            plt.figure(figsize=(10, 6))
            gene_data = gene_expression_df.loc[genes_in_cluster].T  # Transpose for better plotting
            gene_data.plot(ax=plt.gca(), title=f'Gene Expression for {current_id}')
            plt.ylabel('Expression Level')
            plt.xlabel('Samples')

            # Step 6.1: Add all sample names on the x-axis explicitly
            plt.xticks(ticks=range(len(gene_data.index)), labels=gene_data.index, rotation=45, ha='right')  # Display all samples
            plt.tight_layout()

            # Step 6.2: Save the plot with a unique filename
            plt.savefig(os.path.join(cluster_folder, f'{current_id}_genes_expression_cluster_{cluster_id}.png'))
            plt.close()

        # Step 7: Plot metabolite abundance pattern for the current ID
        if metabolites_in_cluster:
            plt.figure(figsize=(10, 6))
            metabolite_data = metabolite_abundance_df.loc[metabolites_in_cluster].T  # Transpose for better plotting
            metabolite_data.plot(ax=plt.gca(), title=f'Metabolite Abundance for {current_id}')
            plt.ylabel('Abundance Level')
            plt.xlabel('Samples')

            # Step 7.1: Add all sample names on the x-axis explicitly
            plt.xticks(ticks=range(len(metabolite_data.index)), labels=metabolite_data.index, rotation=45, ha='right')  # Display all samples
            plt.tight_layout()

            # Step 7.2: Save the plot with a unique filename
            plt.savefig(os.path.join(cluster_folder, f'{current_id}_metabolites_abundance_cluster_{cluster_id}.png'))
            plt.close()
