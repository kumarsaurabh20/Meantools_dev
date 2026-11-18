#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import sys
import argparse
import subprocess as sp
from multiprocessing import Pool
from itertools import combinations
import pandas as pd
import numpy as np
from scipy import stats
from tqdm import tqdm



def calc_MR(rank1, rank2):
    MR = np.sqrt(rank1 * rank2)
    return MR


def get_MR_from_correlations(correlations_df, all_objects):
    # GET RANKS
    ranks_dict = {}
    for cur_object in tqdm(all_objects, total=len(all_objects), desc='Getting ranks'):
        df1 = correlations_df[correlations_df.metabolite == cur_object]
        df1 = df1[['gene', 'correlation']]
        df1.columns = ['other_object', 'correlation']

        df2 = correlations_df[correlations_df.gene == cur_object]
        df2 = df2[['metabolite', 'correlation']]
        df2.columns = ['other_object', 'correlation']

        dff = pd.concat([df1, df2]).sort_values(by='correlation', ascending=False)
        dff.index = np.arange(1, len(dff) + 1)  # this makes index = rank(cur_gene -> other_gene)

        ranks_dict[cur_object] = {}
        for rank, cur_row in dff.iterrows():
            ranks_dict[cur_object][cur_row.other_object] = rank

    # GET MUTUAL RANKS
    edges_df = []
    object_combinations = {n for n in combinations(all_objects, 2)}
    for metabolite, gene in tqdm(object_combinations, total=len(object_combinations), desc='Getting mutual ranks'):
        if ((metabolite in ranks_dict and gene in ranks_dict and gene in ranks_dict[metabolite] and metabolite in ranks_dict[gene])):
            rank1 = ranks_dict[metabolite][gene]
            rank2 = ranks_dict[gene][metabolite]
            edges_df.append([metabolite, gene, calc_MR(rank1, rank2)])
    edges_df = pd.DataFrame.from_records(edges_df, columns=['metabolite', 'gene', 'MR'])

    # OUTPUT
    #fname = os.path.join(Options.output_folder, 'MR.csv')
    #edges_df.to_csv(fname, index=False)
    return edges_df

def get_weight_from_MR(edges_df, edge_weight_cutoff, decay_rate):
    # CONVERT TO EDGE WEIGHT
    print('Estimating edge weights using exponetial decay...')
    edges_df['weight'] = np.exp(-(edges_df.MR - 1) / decay_rate)  # other options: /5 /10 /25 /50 /100
    weights_df = edges_df[edges_df.weight > edge_weight_cutoff].reset_index(drop=True) # default edge_weight_cutoff is 0.01
    
    return weights_df

def validate_clusterone(clusterone):
    if not os.path.isfile(clusterone):
        raise RuntimeError(f"ClusterOne JAR not found: '{clusterone}'.")
    try:
        result = sp.run(['java', '-version'], stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError("Java was not found. Install a Java runtime in the active environment "
                           "and verify that 'java -version' succeeds.") from exc
    if result.returncode:
        raise RuntimeError("Java runtime is unavailable. Install Java and verify 'java -version'. "
                           + (result.stderr or result.stdout).strip())


def run_clusterone(outfile, clusterone, weights_file):
    validate_clusterone(clusterone)
    try:
        with open(outfile, 'w') as f:
            sp.run(['java', '-jar', clusterone, weights_file, '-F', 'csv'],
                   stdout=f, stderr=sp.PIPE, text=True, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError("Java or a ClusterOne input file was not found: " + str(exc)) from exc
    except sp.CalledProcessError as exc:
        diagnostic = (exc.stderr or '').strip()
        hint = (" Check for whitespace in node identifiers; ClusterOne expects two node IDs "
                "followed by a numeric weight." if 'NumberFormatException' in diagnostic else '')
        raise RuntimeError(f"ClusterOne failed (exit {exc.returncode}).{hint}\n{diagnostic}") from exc
    if os.path.getsize(outfile) == 0:
        raise RuntimeError(f"ClusterOne produced empty output: '{outfile}'. Check the input network.")

