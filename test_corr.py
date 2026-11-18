#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def calculate_correlation(df1, df2):
    correlations = []
    for _, row1 in df1.iterrows():
        row_correlations = []
        for _, row2 in df2.iterrows():
            correlation = np.corrcoef(row1.sort_values(), row2.sort_values())[0, 1]
            row_correlations.append(correlation)
        correlations.append(row_correlations)
    return correlations

def save_correlation_images(correlations):
    for i, row_correlations in enumerate(correlations):
        plt.figure()
        plt.plot(range(len(row_correlations)), row_correlations)
        plt.xlabel('Row Index')
        plt.ylabel('Correlation')
        plt.title(f'Correlation for Row {i}')
        plt.savefig(f'correlation_row_{i}.png')
        plt.close()

# Example usage
df1 = pd.read_csv("simul_data_tr.csv", index_col=0)
print(df1)
df2 = pd.read_csv("simul_data_ms.csv", index_col=0)
print(df2)

df1_labels = set(df1.columns)
df2_labels = set(df2.columns)
common_labels = sorted(list(df1_labels.intersection(df2_labels)))
df1 = df1[common_labels]
df2 = df2[common_labels]

correlations = calculate_correlation(df1, df2)
save_correlation_images(correlations)