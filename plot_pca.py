import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Load the expression matrix and metadata file
expression_file = '/Users/singh018/Documents/Meantools_v1/Jeon_tomato/rnaseq.rpkm.csv'
metadata_file = 'TR_pheno.csv'

# Load the expression matrix and metadata
expression_df = pd.read_csv(expression_file, index_col=0)
metadata_df = pd.read_csv(metadata_file)

# Step 1: Standardize the expression data (genes as rows, samples as columns)
scaler = StandardScaler()
expression_scaled = scaler.fit_transform(expression_df.T)  # Transpose to have samples as rows

# Step 2: Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(expression_scaled)

# Step 3: Create a PCA DataFrame
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'], index=expression_df.columns)

# Merge PCA results with metadata
pca_with_metadata = pca_df.merge(metadata_df, left_index=True, right_on='sample')

# Step 4: Map colors for conditions and shapes for time points
# Automatically generate colors and shapes
unique_conditions = pca_with_metadata['condition'].unique()
unique_timepoints = pca_with_metadata['time'].unique()

# Use seaborn's color palette for distinct colors and different marker shapes
colors = sns.color_palette('husl', len(unique_conditions))
color_map = dict(zip(unique_conditions, colors))

markers = ['o', 's', '^', 'D', 'v', 'P', '*']  # Marker shapes for different time points
marker_map = dict(zip(unique_timepoints, markers))

# Step 5: Plot the PCA results
plt.figure(figsize=(10, 8))

# Loop through each unique combination of condition and time point
for condition in unique_conditions:
    for time in unique_timepoints:
        subset = pca_with_metadata[(pca_with_metadata['condition'] == condition) &
                                   (pca_with_metadata['time'] == time)]
        
        plt.scatter(subset['PC1'], subset['PC2'], 
                    color=color_map[condition], 
                    marker=marker_map[time], 
                    label=f'{condition} - {time}h', 
                    alpha=0.7)

# Customize the plot
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance explained)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance explained)')
plt.title('PCA of Gene Expression (Colored by Condition, Shaped by Time)')
plt.legend(loc='best', bbox_to_anchor=(1.05, 1), fontsize=10)
plt.tight_layout()
plt.show()



## ONLY TIME

'''
# Step 4: Automatically generate colors for time points instead of conditions
# Generate a color palette for the unique time points
time_colors = sns.color_palette('husl', len(unique_timepoints))
time_color_map = dict(zip(unique_timepoints, time_colors))

# Plot the PCA with colors representing time points
plt.figure(figsize=(10, 8))

# Loop through each unique time point and plot with respective colors
for time in unique_timepoints:
    subset = pca_with_metadata[pca_with_metadata['time'] == time]
    
    plt.scatter(subset['PC1'], subset['PC2'], 
                color=time_color_map[time], 
                label=f'Time {time}h', 
                alpha=0.7)

# Customize the plot
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance explained)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance explained)')
plt.title('PCA of Gene Expression (Colored by Time)')
plt.legend(loc='best', bbox_to_anchor=(1.05, 1), fontsize=10)
plt.tight_layout()
plt.show()

'''


