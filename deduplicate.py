import pandas as pd

# Sample DataFrame
#data = {'id': ['mc1'],
#        'metabolites': ['ABC, BVG, ABC, HGH, HHG, HGH, JJU, JJY, JJT'],
#        'genes': ['djd, djdjd, jdndj, djd, jhy, jhy']}

#df = pd.read_csv("Jeon_bac_N_10_cluster_formatted.csv", delimiter=',', quotechar='"')

m = pd.read_csv("metabolites.csv")
g = pd.read_csv("genes.csv")

print(m)
print(g)

m = m.drop_duplicates(subset=['metabolites'], keep='first')
g = g.drop_duplicates(subset=['genes'], keep='first')

print(m)
print(g)

m_string = ', '.join(m['metabolites'])
g_string = ', '.join(g['genes'])

final_df = pd.DataFrame({'id': 'mc1', 'metabolites': [m_string], 'genes': [g_string]})
print(final_df)

#final_df.to_csv("merged_clusters_metoverlapped_deduplicated.csv", index=False)

# Aggregate the values back into comma-separated strings if needed
#df = df.groupby('id').agg({
#    'metabolites': lambda x: ', '.join(x),
#    'genes': lambda x: ', '.join(x)
#}).reset_index()

#print(df)
