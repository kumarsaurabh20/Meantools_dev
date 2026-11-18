import os
import sys
import re
import glob
import argparse
#
import pandas as pd
from tqdm import tqdm
from scipy.stats import zscore
import sqlite3
import seaborn as sns
from matplotlib.colors import to_hex
#
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
#
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.cluster import SpectralClustering, Birch, AffinityPropagation
from sklearn.cluster import MeanShift, estimate_bandwidth
#
import gizmos


#python test_time_split_FC.py \
#-ft /Users/singh018/Documents/Meantools_v1/Jeon_tomato/new_MS_data/new_MS_data_combined.csv 
#-qm /Users/singh018/Documents/Meantools_v1/Jeon_tomato/rnaseq.rpkm.csv 
#-dr 25 
#-dn /Users/singh018/Documents/Meantools_v1/Jeon_tomato/Jeon_bac.sqlite


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-ft','--feature_table', default='', action='store', required=True, help='Metabolomics-based feature table!')
	parser.add_argument('-qm','--quantitation_matrix', default='', action='store', required=True, help='Normalized expression table from the RNAseq data!')
	parser.add_argument('-dr', '--decay_rate', default=25, type=int, required=False, help='Deacy rate of FCs. Default: 25')
	parser.add_argument('-g', '--groups', default='treatment', required=False, choices=['treatment', 'timepoint'], help='Group samples based either on timepoint or treatment!')
	parser.add_argument('-m','--metadata', default='', action='store', required=False, help='metadata file representing your experimental design!')
	parser.add_argument('-o', '--outfile', default=False, required=False, action='store', help='Name of the matrix file showing association strength (WARNING:: file size is dependent on number of FCs)!')
	parser.add_argument('-dn', '--sqlite_db_name', default='', action='store', required=True, help='Provide a name for the database')
	return parser.parse_args()



def combine_and_plot_pca(fingerprints_by_group, label_name, output_prefix):
    """
    combine fingerprints from all groups and generate a combined PCA plot.

    Parameters:
    - fingerprints_by_group: Dictionary of fingerprints DataFrames for each group
    - label_name: The label for the grouping (e.g., 'Time' or 'Treatment')
    - output_prefix: Prefix for output files
    """
    combined_df = pd.DataFrame()

    for group, fingerprints_df in fingerprints_by_group.items():
        fingerprints_df['Group'] = group  #add group label to each DataFrame
        combined_df = pd.concat([combined_df, fingerprints_df], ignore_index=True)

    #perform PCA on combined data
    pca_df, explained_variance = pca_on_metafingerprints(combined_df)

    #add group labels to PCA DataFrame
    pca_df['Group'] = combined_df['Group']

    #save combined PCA results
    combined_pca_output = f"{output_prefix}_combined_pca.csv"
    pca_df.to_csv(combined_pca_output, index=False)
    print(f"Saved combined PCA results to {combined_pca_output}")

    #plot combined PCA
    plt.figure(figsize=(10, 8))
    unique_groups = pca_df['Group'].unique()
    for group in unique_groups:
        group_data = pca_df[pca_df['Group'] == group]
        plt.scatter(group_data['PC1'], group_data['PC2'], label=group, alpha=0.7)

    plt.title(f'Combined PCA Plot ({label_name})')
    plt.xlabel(f'Principal Component 1 ({explained_variance[0]:.2f}% variance explained)')
    plt.ylabel(f'Principal Component 2 ({explained_variance[1]:.2f}% variance explained)')
    plt.legend(title=label_name)
    plt.grid(True)
    combined_pca_plot = f"{output_prefix}_combined_pca_plot.png"
    plt.savefig(combined_pca_plot, dpi=300)
    print(f"Saved combined PCA plot to {combined_pca_plot}")
    plt.show()


def group_samples(metadata_df, transcripts_df, metabolites_df, group_by='time', group_values=None):
    
    """
    group samples based on time points or treatments and filter corresponding data.

    Parameters:
    - metadata_df: DataFrame containing metadata
    - transcripts_df: Transcriptomics DataFrame
    - metabolites_df: Metabolomics DataFrame
    - group_by: Grouping criteria ('time' or 'condition')
    - group_values: List of values to group by (e.g., [12, 24, 48] for time or ['M', 'F'] for treatments)

    Returns:
    - a dictionary with keys as groups (e.g., time points or treatment groups) and filtered DataFrames
    """
    grouped_data = {}
    for group in group_values:
        if group_by == 'time':
            subset = metadata_df[metadata_df['time'] == group]
        elif group_by == 'condition':
            subset = metadata_df[metadata_df['condition'] == group]
        else:
            raise ValueError("Invalid group_by value. Use 'time' or 'condition'.")
        
        selected_samples = subset['sample'].values
        transcripts_subset = transcripts_df[['feature'] + [col for col in transcripts_df.columns if col in selected_samples]]
        metabolites_subset = metabolites_df[['feature'] + [col for col in metabolites_df.columns if col in selected_samples]]
        grouped_data[group] = (transcripts_subset, metabolites_subset)
    
    return grouped_data


def calculate_fingerprints_by_group(grouped_data, clusters_frame):
    """
    calculate fingerprints for each group.

    parameters:
    - grouped_data: Dictionary of grouped transcriptomics and metabolomics data
    - clusters_frame: DataFrame with cluster information

    Returns:
    - dictionary of fingerprints for each group
    """
    fingerprints_by_group = {}

    for group, (transcripts_df, metabolite_df) in grouped_data.items():
        results_rows = []
        for _, row in tqdm(clusters_frame.iterrows(), total=clusters_frame.shape[0], desc=f'Processing Group {group}'):
            results_rows.append(calc_metafingerprints(row, transcripts_df, metabolite_df))
        
        fingerprints_df = pd.concat(results_rows, axis=1).T
        fingerprints_by_group[group] = fingerprints_df
    
    return fingerprints_by_group


def perform_clustering_and_extract_ids(df, x, y, n_clusters, method='', **kwargs):
    
	"""
    apply a clustering method (KMeans, DBSCAN, Agglomerative) on PCA-reduced data and add cluster labels to the DataFrame.
    
    Parameters:
    - pca_df: DataFrame containing PCA components and 'ID' column
    - method: The clustering method to use ('kmeans', 'dbscan', 'agglomerative')
    - x, y: The principal components to cluster on (e.g., PC1 and PC2)
    - kwargs: Additional parameters specific to each clustering method
    
    Returns:
    - Updated DataFrame with an added 'Cluster' column
    """
	if method == 'kmeans':
		n_clusters = kwargs.get('n_clusters', 2)  # Default to 3 clusters if not provided
		model = KMeans(n_clusters=n_clusters, random_state=42)

	elif method == 'dbscan':
		eps = kwargs.get('eps', 0.5)  # Default epsilon value if not provided
		min_samples = kwargs.get('min_samples', 5)  # Default minimum samples if not provided
		model = DBSCAN(eps=eps, min_samples=min_samples)

	elif method == 'agglomerative':
		n_clusters = kwargs.get('n_clusters', 3)  # Default to 3 clusters if not provided
		model = AgglomerativeClustering(n_clusters=n_clusters)
	
	elif method == 'gmm':
		n_components = kwargs.get('n_components', 3)
		model = GaussianMixture(n_components=n_components, random_state=42)

	elif method == 'spectral':
		n_clusters = kwargs.get('n_clusters', 3)
		model = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors', random_state=42)

	elif method == 'meanshift':
		bandwidth = estimate_bandwidth(df[[f'PC{x}', f'PC{y}']])
		model = MeanShift(bandwidth=bandwidth)

	elif method == 'birch':
		n_clusters = kwargs.get('n_clusters', 3)
		model = Birch(n_clusters=n_clusters)
		pca_data = np.ascontiguousarray(df[[f'PC{x}', f'PC{y}']].values)
		df['Cluster'] = model.fit_predict(pca_data)
		return df

	elif method == 'affinity':
		damping = kwargs.get('damping', 0.9)
		preference = kwargs.get('preference', -50)
		model = AffinityPropagation(damping=damping, preference=preference)

	else:
		raise ValueError(f"Unknown method '{method}'. Supported methods: 'kmeans', 'dbscan', 'agglomerative'.")

	#fit the model to the PCA-reduced data (PCx and PCy)
	df['Cluster'] = model.fit_predict(df[[f'PC{x}', f'PC{y}']])

	# Step 1: Assign colors to clusters
	unique_clusters = df['Cluster'].unique()
	color_palette = sns.color_palette('hsv', len(unique_clusters))  # Assign distinct colors
	cluster_colors = dict(zip(unique_clusters, color_palette))

	# Step 2: Add a column to store the cluster color (in hex format)
	df['Cluster_Color'] = df['Cluster'].map(lambda c: to_hex(cluster_colors[c]))


	return df



def calculate_sigma(df):
	# Calculate the singular value decomposition (SVD) for the DataFrame 'df'
	# Here, I'm using 'Column1' and 'Column2' as the columns for the SVD calculation
	# Replace them with your actual column names
	numeric_cols = df.select_dtypes(include=[np.number]).columns
	X = df[numeric_cols].values
	U, s, VT = np.linalg.svd(X, full_matrices=False)
	#sigma = np.diag(s) ## use this if you need a matrix instead of a 1D vector
	sigma = VT[0]	## VERY IMPORTANT::First right singular vector
	return sigma



def fill_extracted_rows(row, transcripts_df, metabolite_df):

	extracted_rows_combined = pd.DataFrame()
	cluster_id = row.iloc[0]
	elements_str = row.iloc[1]
	elements_list = elements_str.split()

	for element in elements_list:
		if element in transcripts_df['feature'].values:
			extracted_row = transcripts_df[transcripts_df['feature'] == element]
			extracted_rows_combined = pd.concat([extracted_rows_combined, extracted_row], ignore_index=True)
		elif element in metabolite_df['feature'].values:
			extracted_row = metabolite_df[metabolite_df['feature'] == element]
			extracted_rows_combined = pd.concat([extracted_rows_combined, extracted_row], ignore_index=True)
		else:
			sys.exit("element is not found in either data frame!")

	return cluster_id, extracted_rows_combined



def calc_metafingerprints(row, transcripts_df, metabolite_df):

	cluster_id, extracted_rows_combined = fill_extracted_rows(row, transcripts_df, metabolite_df)

	#return the right singular matrix by decomposing a matrix using svd.
	sigma = calculate_sigma(extracted_rows_combined)

	row['metafingerprints'] = sigma
	return row

def plot_pca(df, explained_variance, x, y, ids_to_highlight=None, cluster_id=False):

	plt.figure(figsize=(8, 6))

	if cluster_id:
		scatter = plt.scatter(df[f'PC{x}'], df[f'PC{y}'], c=df['Cluster'], cmap='viridis', alpha=0.6)
		plt.colorbar(scatter, label='Cluster')
	else:
		plt.scatter(df[f'PC{x}'], df[f'PC{y}'], alpha=0.6)
	
	if ids_to_highlight is not None:
		for i, row in df.iterrows():
			if row['ID'] in ids_to_highlight:
				plt.annotate(row['ID'], (row[f'PC{x}'], row[f'PC{y}']),textcoords="offset points", xytext=(5,-5), ha='right', fontsize=9, color='red')


	plt.title('PCA of Metafingerprints')
	plt.xlabel(f'Principal Component {x} ({explained_variance[x-1]:.2f}% variance explained)')
	plt.ylabel(f'Principal Component {y} ({explained_variance[y-1]:.2f}% variance explained)')
	plt.show()


def pca_on_metafingerprints(df):

	metafingerprints = np.array(df['metafingerprints'].tolist())
	scaler = StandardScaler()
	metafingerprints_scaled = scaler.fit_transform(metafingerprints)
	pca = PCA(n_components=5)
	
	principal_components = pca.fit_transform(metafingerprints_scaled)
	explained_variance = pca.explained_variance_ratio_ * 100
	pca_df = pd.DataFrame(data=principal_components, columns=[f'PC{i+1}' for i in range(principal_components.shape[1])])
	pca_df['ID'] = df['ID'].values
	plot_pca(pca_df, explained_variance, 1, 2)
	#plot_pca(pca_df, explained_variance, 3, 4)
	#plot_pca(pca_df, explained_variance, 4, 5)

	return pca_df, explained_variance


def main(Options):

	#raw expression and features file.
	#used for computing fingerprints
	transcripts = pd.read_csv(Options.quantitation_matrix)
	transcripts = transcripts.rename(columns={'gene': 'feature'})
	metabolites = pd.read_csv(Options.feature_table)
	metabolites = metabolites.rename(columns={'metabolite': 'feature'})

	#keep common labels between the two datasets
	#it is important to keep the gene and metabolite column names (from the main omics dataset) 
	#as 'feature'
	transcripts_labels = set(transcripts.columns)
	metabolites_labels = set(metabolites.columns)
	common_labels = sorted(list(transcripts_labels.intersection(metabolites_labels)))
	transcripts_df = transcripts[common_labels]
	metabolite_df = metabolites[common_labels]
	
	#scale the transcriptomics and metabolomics datasets by using z-scores 
	#especially for heatmaps
	t_numeric_columns = transcripts_df.copy()
	t_numeric_columns[t_numeric_columns.select_dtypes(['float64', 'int64']).columns] = t_numeric_columns.select_dtypes(['float64', 'int64']).apply(zscore)
	transcripts_df = t_numeric_columns
	#
	m_numeric_columns = metabolite_df.copy()
	m_numeric_columns[m_numeric_columns.select_dtypes(['float64', 'int64']).columns] = m_numeric_columns.select_dtypes(['float64', 'int64']).apply(zscore)
	metabolite_df = m_numeric_columns

	metadata = pd.read_csv(Options.metadata)

    #parse group selection flags
	if Options.groups == "timepoint":
		grouped_data = group_samples(metadata, transcripts, metabolites, group_by='time', group_values=['12','24','48'])
	elif Options.groups == "treatment":
		grouped_data = group_samples(metadata, transcripts, metabolites, group_by='condition', group_values=['M', 'F', 'B'])
	else:
		raise ValueError("You must specify --groups!")


	#read the cluster files with multiple DR
	all_files = gizmos.import_from_sql(Options.sqlite_db_name, sqlite_tablename="", df_columns=[], conditions={}, structures = False, clone = True)
	
	if not all_files:
		raise ValueError(f"No data was imported from {Options.sqlite_db_name}.")

	#adding the DR info to know which cluster is coming from which file
	#DR is deacy rate used to calculate clusters
	for table_name, df in all_files.items():
		match = re.search(r'(DR_\d+)', table_name)
		df['Source'] = match.group(1) if match else 'unknown'
		#if match:
		#	source_value = match.group(1)
		#	df['Source'] = source_value
		#else:

	#concatenate all DataFrames into a single DataFrame
	frame = pd.concat(all_files.values(), axis=0, ignore_index=True)

	list_of_counts = []
	for i in range(1, frame.shape[0] + 1):
		list_of_counts.append(f'cluster{i:03d}')
	frame['cluster_id'] = list_of_counts


	#merge cluster ID and source info to create a new ID
	#frame["ID"] = frame["cluster_id"].str.cat(frame["Source"], sep = "_")
	frame["ID"] = frame["cluster_id"] +"_"+ frame["Source"] + "_" + frame["Cluster"].astype(str)


	clusters_frame = frame[['ID', 'Members']]
	clusters_frame.to_csv("Ara_clusters_members.csv", index=False)

	filtered_df = clusters_frame[clusters_frame['ID'].str.contains('DR_25')]


	clusters_with_fingerprints_df = pd.DataFrame()
	results_rows = []


	#calculate fingerprints for each group
	fingerprints_by_group = calculate_fingerprints_by_group(grouped_data, clusters_frame)
	
	#save results for each group
	for group, fingerprints_df in fingerprints_by_group.items():
		output_file = f"fingerprints_group_{group}.csv"
		fingerprints_df.to_csv(output_file, index=False)
		print(f"Saved fingerprints for group '{group}' to {output_file}")

		#perform PCA
		pca_df, explained_variance = pca_on_metafingerprints(fingerprints_df)

		#cluster PCA results
		pca_df = perform_clustering_and_extract_ids(pca_df, 1, 2, 3, method='kmeans')
		pca_output_file = f"pca_group_{group}_labelled.csv"
		pca_df.to_csv(pca_output_file, index=False)
		print(f"Saved PCA results for group '{group}' to {pca_output_file}")

		#plot PCA
		plot_pca(pca_df, explained_variance, 1, 2, ids_to_highlight=None, cluster_id=True)



'''
	for idx, row in tqdm(filtered_df.iterrows(), total=filtered_df.shape[0], desc='Processing FCs'):
		results_rows.append(calc_metafingerprints(row, transcripts_df, metabolite_df))
	clusters_with_fingerprints_df = pd.concat(results_rows, axis=1).T

	#print(clusters_with_fingerprints_df)

	pca_df, explained_variance = pca_on_metafingerprints(clusters_with_fingerprints_df)

	pca_df = perform_clustering_and_extract_ids(pca_df, 1, 2, 3, method='kmeans')
	pca_df.to_csv("tomato_clusters_DR_25_pca_labelled.csv", index=False)

	plot_pca(pca_df, explained_variance, 1, 2, ids_to_highlight=None, cluster_id=True)
'''



if __name__ == "__main__":

	Options = get_args()
	main(Options)