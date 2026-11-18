import pandas as pd


def split_and_match(row, metabolite_df, column):

	'''
	The function splits elements and return matching/unmatching
	Splits genes and metabolites into separate columns

	:param1 row: row of the dataframe
	:param2 metabolite_df: because it checks metabolic features within a merged_cluster, unmatched elements automatically assigned to genes column
	
	:return : returns two combined series which can be integrated into a data frame
	'''

	merged_members = row[column].split(',')
	matching = [metabolite for metabolite in merged_members if metabolite in metabolite_df['metabolite'].to_list()]
	unmatching = [metabolite for metabolite in merged_members if metabolite not in matching]
	return pd.Series([matching, unmatching])


list_df = pd.read_csv("Jeon_bac_cluster_N_10.csv", quotechar='"')
print(list_df)
metabolite_df = pd.read_csv("MS.combined.metabolome.csv")
print(metabolite_df)
cluster_df = pd.DataFrame()
cluster_df[['metabolites', 'genes']] = list_df.apply(split_and_match, args=(metabolite_df, 'Members'), axis=1)
print(cluster_df)
# Convert the list in the 'Cluster' column to a string
cluster_df['metabolites'] = cluster_df['metabolites'].apply(lambda x: ', '.join(x))
cluster_df['genes'] = cluster_df['genes'].apply(lambda x: ', '.join(x))
cluster_df.to_csv('Jeon_bac_N_10_cluster_formatted.csv', index=False)

