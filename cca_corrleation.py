import pandas as pd
import numpy as np
#from scipy.stats import canonical_correlation_analysis
from sklearn.cross_decomposition import CanonicalCorrelationAnalysis

# Create dummy data for the first dataset
data1 = {
    'gene': ['CRO-01G000010', 'CRO-01G000020'],
    'idioblast': [0, 0.047381830244514],
    'epidermis': [0.0207128964707296, 0.183657899150503]
}
df1 = pd.DataFrame(data1)

# Create dummy data for the second dataset
data2 = {
    'metabolite': ['CR0025', 'CR0024'],
    'idioblast': [0.0335479223991636, 0.00898131204583506],
    'epidermis': [0.0989849959455104, 0.0205645402698884]
}
df2 = pd.DataFrame(data2)

# Extract the "idioblast" and "epidermis" columns from both datasets
X = df1[['idioblast', 'epidermis']]
Y = df2[['idioblast', 'epidermis']]

# Convert the DataFrames to numpy arrays
X_array = X.to_numpy()
Y_array = Y.to_numpy()

# Perform Canonical Correlation Analysis
cca = CanonicalCorrelationAnalysis(n_components=2)
cca.fit(X, Y)


# Get canonical correlations and coefficients
canonical_corr = cca.canonical_correlations_
coefficients_x = cca.x_weights_
coefficients_y = cca.y_weights_

# Display canonical correlations and coefficients
print("Canonical Correlations:", canonical_corr)
print("Canonical Coefficients (X):", coefficients_x)
print("Canonical Coefficients (Y):", coefficients_y)








