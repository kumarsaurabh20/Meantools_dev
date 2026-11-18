import pandas as pd

# Load the first CSV file (root_7dpi_MGH_MC1.csv)
root_df = pd.read_csv('root_7dpi_MGH_MC1.csv')
print(root_df)

# Load the second CSV file (Ara.gene.annotations.csv)
annotations_df = pd.read_csv('/Users/singh018/Downloads/NWO_groot_paired_omics_2024/annotation/Ara.genes.pfam.csv')
print(annotations_df)

genes_list = root_df['genes'][0].split(",")
genes_df = pd.DataFrame(genes_list, columns=['gene'])
print(genes_df)

# Merge on the 'genes' column to get descriptions
merged_df = pd.merge(genes_df, annotations_df, left_on='gene', right_on='transcript_id', how='left')
merged_df = merged_df.drop_duplicates()

print(merged_df)
# Save the output to a new CSV file

merged_df.to_csv('MC1_genes.pfam.csv', index=False)

print("Merged file saved as 'merged_output.csv'.")
