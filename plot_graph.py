from random import sample
from pyvis.network import Network
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import sqlite3


def import_from_sql(sqlite_db_name="", sqlite_tablename="", conditions={}):

    if os.path.isfile(sqlite_db_name):
        conn = sqlite3.connect(sqlite_db_name)
    else:
        raise FileNotFoundError("The {} file was not found! Did you run the correlation step first?".format(sqlite_db_name))


    # Connect to the SQLite database
    conn = sqlite3.connect(sqlite_db_name)
    
    df_columns = ['metabolite', 'gene', 'MR', 'weight']
        
    query = f"SELECT {', '.join(df_columns)} FROM {sqlite_tablename} WHERE {' AND '.join([f'{key} IN {tuple(conditions[key])}' for key in conditions])};"

    # Execute the query and fetch the results
    cursor = conn.cursor()
    cursor.execute(query)
    results = cursor.fetchall()

    # Display or further process the extracted data
    result_df = pd.DataFrame(results, columns=df_columns)

    conn.close()

    return result_df


### Start here

meregd_clusters_df = pd.read_csv("2HoTRe_cluster256.csv", index_col='id', delimiter=',', quotechar='"')

meregd_clusters_df['genes'] = meregd_clusters_df['genes'].apply(lambda x: ','.join(set(x.split(','))))
meregd_clusters_df['metabolites'] = meregd_clusters_df['metabolites'].apply(lambda x: ','.join(set(x.split(','))))
meregd_clusters_df = meregd_clusters_df.drop_duplicates(subset=['genes', 'metabolites'], keep="first")

print(meregd_clusters_df)

temp = meregd_clusters_df.to_dict(orient='list')

meregd_clusters_dict = {}
for key, values in temp.items():    
    new_key = key[:-1]  # Remove the trailing 's' from the key
    new_values = [item.strip() for value in values for item in value.split(',')]  # Split and flatten the values
    meregd_clusters_dict[new_key] = new_values
    
correlation_df = import_from_sql(sqlite_db_name="roots_7dpi_mgh_met_up.sqlite", sqlite_tablename="roots_7dpi_met_UP_MR_weights_DR_25", conditions=meregd_clusters_dict)
#correlation_df = import_from_sql(sqlite_db_name="/Users/singh018/Downloads/NWO_groot_paired_omics_2024/NWO_groot_coexpression.sqlite", sqlite_tablename="nwo_combined_MR_weights_DR_25", conditions=meregd_clusters_dict)

df = correlation_df[['metabolite','gene', 'weight']]
df = df.dropna()
print(df)

#m = {'M258T1128':1, 'M274T1158':2, 'M277T882':2, 'M600T895':3, 'M274T851':2, 'M276T882':4}
#pal
#m={'M468T649':1, 'M469T650':1, 'M470T650':1, 'M284T651':1, 'M243T650':1, 'M446T648':1}
#data = {'metabolite':['Ar329R754', 'Ar327R754', 'Ar679R444', 'Ar282R1153', 'Ar445R1094', 'Ar403R1059', 'Ar637R1131'], 'class_group':[1,1,1,1,1,1,1]}
#THT
#data = {'metabolite':['M316T772'], 'class_group':[1]}
#falca
#data = {'metabolite': ['M258T1128','M274T1158','M277T882','M600T895','M274T851','M276T882'], 'class_group':[1,1,1,1,1,1]}
#m_df = pd.DataFrame(data)
#print(m_df)
#merged_df = pd.merge(df, m_df, on="metabolite", how="left")
#print(merged_df)

df['class_group'] = df['metabolite'].apply(lambda x: 1)
print(df)
#df.to_csv("combined_cluster_edge_weights.csv", index=False)


#PAL
#1: Isoflavanoids (Superclass: Phenylpropanoids and polyketides)
#2: 2-arylbenzofuran flavanoids (Superclass: Phenylpropanoids and polyketides)
#3: Aurone flavanoids ((Superclass: Phenylpropanoids and polyketides))

#Falcarindiol
#1: stilbenes (Superclass: phenylpropanoids and polyketides)
#2: prenol lipids (Superclass = Lipid-like molecules; subclass = sesquiterpenoids)
#3: Steroid like drivatives (Superclass = Lipid-like derivatives; subclass = pregane steroids or steroid esters)
#4: Carboxyl acid derivatives (Superclass = Organic acids; subclass = tetra carboxylic acids
#Organooxygen compounds (Superclass = Organic oxygen compounds; subclass = Alcohols and polyols)

g_from_data = Network(height='95%',width='90%',
              bgcolor='white',font_color="black",
              heading="New feature network based on MR N25")

destinations = list(df['metabolite'])
sources = list(df['gene'])
weights = list(df['weight'])
classes = list(df['class_group'])

edge_data = zip(sources, destinations, weights, classes)

for e in edge_data:
                src = e[0]
                dst = e[1]
                w = e[2]
                c = e[3]

                g_from_data.add_node(src, src, title=src, color='#73c375', size = 30)
                g_from_data.add_node(dst, dst, title=dst, group=c, color='#7777f8')
                g_from_data.add_edge(src, dst, value=w, color='#d3d3d3')

neighbor_map = g_from_data.get_adj_list()

for node in g_from_data.nodes:
                node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
                node["value"] = len(neighbor_map[node["id"]])

g_from_data.show_buttons(filter_ = ['nodes', 'edges', 'physics'])
g_from_data.show('2HoTRe_cluster256_N25.html')
#net = Network(notebook=True)
#G = nx.from_pandas_edgelist(file_pd, source='metabolite', target='gene', edge_attr='correlation_coefficient')



#g_from_data = nx.Graph()

#for i in range(len(sources)):
#  try:
#    g_from_data.add_node(sources[i],label=sources[i], color='green')
#  except:
#    pass

#for i in range(len(destinations)):
#  try:
#    g_from_data.add_node(destinations[i],label=destinations[i])
# except:
#   pass

#for i in range(len(second_destinations)):
#  try:
#    g_from_data.add_node(second_destinations[i], label=second_destinations[i], title=smiles[i], color='yellow')  
#  except:
#    pass
 
#for (i,j) in zip(sources, destinations):
#    
#    gene, metabolite = (i,j)
#    try:
#     g_from_data.add_edge(i,j, weight=float(df[(df['metabolite'] == metabolite) & (df['gene'] == gene)]['weight']))
#   except:
#      pass

#for (i,j) in zip(destinations, second_destinations):
#  try:
#    g_from_data.add_edge(i,j)
#  except:
#    pass

#neighbor_map=got_net.get_adj_list()

# add neighbor data to node hover data
#for node in g_from_data.nodes:
#    node['title'] += ' Neighbors:<br>' + '<br>'.join(neighbor_map[node['id']])
#    node['value'] = len(neighbor_map[node['id']])
 
#g_from_data.show_buttons(['physics'])
#pos = nx.spring_layout(g_from_data, seed=0)
#nx.draw_networkx(g_from_data, pos)
#for edge in g_from_data.edges(data=True):
#    nx.draw_networkx_edges(g_from_data, pos, edgelist=[edge], width=edge[2])


#display(HTML('A_Complete_Networkx_Graph_From_DataFrame.html'))

#plt.figure(figsize=(8, 6))
#nx.draw(g_from_data, pos=pos, node_size=10)
#plt.show()


