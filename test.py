from __future__ import division, print_function
import os
import sys
import glob
import argparse
import re
import pandas as pd
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
import svgutils.transform as sg
from svgelements import *

import gizmos
import svg




def get_structure_network_edge_info(cur_edge, df):
    substrate_mask = cur_edge.predicted_substrate_id == df['predicted_substrate_id']
    product_mask = cur_edge.predicted_product_id == df['predicted_product_id']
    # CORRELATED GENES
    cur_edge['genes_corr_substrate'] = len(
        df[substrate_mask & product_mask].dropna(subset=['correlation_substrate']))
    cur_edge['genes_corr_product'] = len(
        df[substrate_mask & product_mask].dropna(subset=['correlation_product']))
    cur_edge['genes_corr_both'] = len(
        df[substrate_mask & product_mask].dropna(subset=['correlation_substrate', 'correlation_product']))

    return cur_edge


def get_structure_network(reactions_df):
    structures_network_df = reactions_df[['predicted_substrate_id', 'predicted_product_id']].drop_duplicates().copy()
    if 'correlation_substrate' in reactions_df:
        structures_network_df = structures_network_df.apply(get_structure_network_edge_info, df=reactions_df, axis=1)
    else:
        structures_network_df['genes_corr_substrate'] = 0
        structures_network_df['genes_corr_product'] = 0
        structures_network_df['genes_corr_both'] = 0
    # ++ genes_corr_substrate, genes_corr_product, genes_corr_both

    structures_network = gizmos.df_to_graph(structures_network_df, use_pd_graph=True)
    # print(structures_network_df)
    #    predicted_substrate_id predicted_product_id genes_corr_substrate genes_corr_product genes_corr_both
    #root
    #NP_27475 NP_27475 VM_3137934 5 5 5
    #NP_27475 NP_27475 VM_7461151 2 2 2

    # print(structures_network)
    #DiGraph with 26 nodes and 25 edges
    return structures_network, structures_network_df


def get_structure_network_attributes(root_id, structures_network, structures_df, reactions_df):
    reaction_distance_df = pd.DataFrame.from_dict(nx.single_source_shortest_path_length(structures_network.to_undirected(), root_id), orient='index', columns=['root_distance'])
    reaction_distance_mod_df = reaction_distance_df.reset_index().rename(columns={'index': 'nodes'})
    list_index = list(reaction_distance_mod_df['nodes'])

    #    root_distance
    #NP_27475   0
    #VM_7461151 1
    #VM_1165316 1
    #VM_4892322 1

    reactions_df_mod = reactions_df.reset_index()
    temp_df = reactions_df_mod[reactions_df_mod.isin(list_index).any(axis=1)]
    substrate = temp_df[['predicted_substrate_id', 'predicted_substrate_smiles']].drop_duplicates('predicted_substrate_id').rename(columns={'predicted_substrate_id': 'predicted_id', 'predicted_substrate_smiles': 'predicted_smiles'})
    substrate['predicted_mm'] = substrate.predicted_smiles.apply(gizmos.get_mm_from_str, is_smarts=True)
    product = temp_df[['predicted_product_id', 'predicted_product_smiles']].drop_duplicates('predicted_product_id').rename(columns={'predicted_product_id': 'predicted_id', 'predicted_product_smiles': 'predicted_smiles'})
    product['predicted_mm'] = product.predicted_smiles.apply(gizmos.get_mm_from_str, is_smarts=True)
    structures_attributes_df = pd.concat([substrate, product], ignore_index=True)
    structures_attributes_dist_df = pd.merge(reaction_distance_mod_df, structures_attributes_df, how="left", left_on='nodes', right_on='predicted_id')
    structures_attributes_dist_df = structures_attributes_dist_df.drop(columns="nodes")
    structures_attributes_dist_df = structures_attributes_dist_df[['predicted_id', 'predicted_smiles', 'predicted_mm', 'root_distance']]
    #print(structures_attributes_dist_df)
    #       predicted_id predicted_smiles predicted_mm root_distance
    #   0   NP_152307 CCCCCCCCC=CCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCCN)... 855.67 0
    #   1   VM_4516106 CCCCCCCCC=CCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCC... 857.69 1
    return structures_attributes_dist_df



gizmos.print_milestone('Loading structures...', verbose=True)
structures_df = pd.read_csv("heraled_HB_filtered/structure_predictions.csv", index_col=0)

if 'predicted_substrate_mm' not in structures_df:
	structures_df['predicted_substrate_mm'] = structures_df.predicted_substrate_smiles.apply(gizmos.get_mm_from_str, is_smarts=False)



gizmos.print_milestone('Loading reactions...', verbose=True)
reactions_df = pd.read_csv("heraled_HB_filtered/reactions.csv", index_col=None, dtype={'rhea_id_reaction': str, 'reaction_id': str})
reactions_df = reactions_df.set_index('root')
            
gizmos.print_milestone('Generating networks...', verbose=True)
            
structures_network, structure_network_df = get_structure_network(reactions_df)

            # SUB NETWORKS
subnetworks = [n for n in nx.connected_components(structures_network.to_undirected())]
subnetworks_dict = dict()
for sub_idx, cur_sub in enumerate(subnetworks):
	for node in cur_sub:
		subnetworks_dict[node] = sub_idx

# PROCESS EACH ROOT NETWORK
gizmos.print_milestone('Processing each network...', verbose=True)
for cur_root in reactions_df.index.unique():

	# GET NETWORK COMPONENT WITH ROOT
	cur_root_structures = subnetworks[subnetworks_dict[cur_root]]
	cur_root_structures_network = structures_network.subgraph(cur_root_structures)


	# GET ROOT-SPECIFIC DATA
	cur_root_structures_df = structures_df[structures_df.index.isin(cur_root_structures)].copy()
    #cur_root_structures_df = structures_df[structures_df['root'].isin(cur_root_structures)]

	subs_mask = structure_network_df.predicted_substrate_id.isin(cur_root_structures)
	prod_mask = structure_network_df.predicted_product_id.isin(cur_root_structures)
	cur_root_structures_network_df = structure_network_df[subs_mask & prod_mask]
	cur_root_structures_attributes_df = get_structure_network_attributes(cur_root, cur_root_structures_network, structures_df, reactions_df)

	root_dag_structures_network = gizmos.get_dag_from_structures_network(cur_root_structures_network, cur_root_structures_attributes_df)
	longest_path = nx.dag_longest_path(root_dag_structures_network)  # list of nodes, ordered

	print(root_dag_structures_network)
	print(longest_path)

