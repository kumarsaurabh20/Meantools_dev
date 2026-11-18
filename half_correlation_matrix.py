from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import gizmos

sns.set_theme(style="white")


#generate a large random dataset
#meregd_clusters_df = pd.read_csv("falcarindiol_cluster.csv", index_col='id', delimiter=',', quotechar='"')
#temp = meregd_clusters_df.to_dict(orient='list')

#meregd_clusters_dict = {}
#for key, values in temp.items():
#    new_key = key[:-1]  # Remove the trailing 's' from the key
#    new_values = [item.strip() for value in values for item in value.split(',')]  # Split and flatten the values
#    meregd_clusters_dict[new_key] = new_values

#corr_columns_for_edge_weight = ['metabolite', 'gene', 'MR', 'weight']
#corr_columns = ['metabolite', 'gene', 'correlation', 'P']
#correlation_df = gizmos.import_from_sql("Jeon_tomato/Jeon_bac.sqlite", "Jeon_bac_correlations", corr_columns, conditions = meregd_clusters_dict, structures = False, clone = False)
#print(correlation_df)
#edge_weight_df = gizmos.import_from_sql("Jeon_tomato/Jeon_bac.sqlite", "Jeon_bac_correlations_MR_weights_DR_10", corr_columns_for_edge_weight, conditions = meregd_clusters_dict, structures = False, clone = False)
#print(edge_weight_df)
#correlation_df = correlation_df[['metabolite', 'gene', 'correlation']]
#edge_weight_df = edge_weight_df[['metabolite', 'gene', 'weight']]
#print(edge_weight_df)
#
correlation_df = pd.read_csv("PAL_correlations_raw.csv")
result = correlation_df.pivot_table(index='metabolite', columns='gene', values='correlation', fill_value=None)
result = result.sort_index(axis=0).sort_index(axis=1)
print(result)

#compute the correlation matrix
#corr = d.corr()

#print(corr)

#generate a mask for the upper triangle
#mask = np.triu(np.ones_like(corr, dtype=bool))

# Set up the matplotlib figure
#f, ax = plt.subplots(figsize=(11, 9))


heatmap_height = len(result) * 0.3
plt.figure(figsize=(15, heatmap_height))

# Generate a custom diverging colormap
#cmap = sns.diverging_palette(230, 20, as_cmap=True)
cmap = sns.color_palette("rocket_r", as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
#sns.heatmap(result, mask=mask, cmap=cmap, vmax=.3, center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})

fig = sns.heatmap(result, annot=True, linewidths=0.5, cmap=cmap, vmin=0, vmax=1) #cbar=False
#plt.xticks(rotation=45)
fig.set(xticklabels=[])
fig.set(xlabel=None)
fig.set(yticklabels=[])
fig.set(ylabel=None)
plt.savefig("Falca_correlations_no_labels.png", bbox_inches='tight', dpi=600)
