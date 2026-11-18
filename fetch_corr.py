#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import argparse
import os
import sys
import numpy as np
import pandas as pd
import sqlite3

import gizmos

sqlite_name = "/Users/singh018/Documents/Meantools_v1/Jeon_tomato/Jeon_bac.sqlite"
tablename = "Jeon_bac_correlations"

meregd_clusters_df = pd.read_csv("/Users/singh018/Documents/Meantools_v1/Jeon_tomato/PAL_4CL_THT_RR/PAL_merged_cluster.csv", index_col='id', delimiter=',', quotechar='"')
temp = meregd_clusters_df.to_dict(orient='list')
    

    # The temp has a dictionary structure but the values are present as string and not as 
    # individual elements.
meregd_clusters_dict = {}
for key, values in temp.items():
	new_key = key[:-1]  # Remove the trailing 's' from the key
	new_values = [item.strip() for value in values for item in value.split(',')]  # Split and flatten the values
	meregd_clusters_dict[new_key] = new_values
    
    #
corr_columns = ['metabolite', 'gene', 'correlation', 'P']
correlation_df = gizmos.import_from_sql(sqlite_name, tablename, corr_columns, conditions = meregd_clusters_dict, structures = False, clone = False)
print(correlation_df)
