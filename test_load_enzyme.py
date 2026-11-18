from __future__ import division, print_function
import os
import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool
from rdkit.Chem import AllChem as Chem

import gizmos


print('Loading coexpression...')
correlation_df = pd.read_csv("Jeon_tomato/Fungal/fungal_correlations.csv", index_col=None)
correlation_df = correlation_df.rename(columns={correlation_df.columns[0]: 'ms_name',
                                                    correlation_df.columns[1]: 'gene',
                                                    correlation_df.columns[2]: 'correlation',
                                                    correlation_df.columns[3]: 'P'})



correlation_df = correlation_df[abs(correlation_df['correlation']) >= 0.7]
correlation_df = correlation_df[abs(correlation_df['P']) <= 0.01]

print("#### This is CORRELATION DF ####")
print(correlation_df)


pfam_dict_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pfams_dict.csv')  # Acc,Name,Desc

    # GENE ANNOTATIONS
    # gene, pfam1;pfam2


print('Loading gene annotations...')
annotations_df = pd.read_csv("Jeon_tomato/tomato.pfams.sol.csv", index_col=None)
#annotations_df = pd.read_csv("gene_annotation.csv", index_col=None)

annotations_df = annotations_df.rename(columns={annotations_df.columns[0]: 'gene', annotations_df.columns[1]: 'enzyme_pfams'})
print(annotations_df)
enzyme_pfams_list = annotations_df.enzyme_pfams.apply(gizmos.pd_to_list, separator=';')
print(enzyme_pfams_list)

    # expand gene annotations so there's one pfam per line (but we keep the "pfam" annotation that have them all)

#lens = [len(str(item).split(",")) for item in enzyme_pfams_list]
lens = [len(item) for item in enzyme_pfams_list]

new_df = pd.DataFrame({'gene': np.repeat(annotations_df.gene, lens), 'pfam_rule': np.concatenate(enzyme_pfams_list)})

print(new_df)

    # Kumar 07/08/2022
    # This step generates a df with 3 columns
    # gene;enzyme_pfams;pfam_rule
    # It makes one gene - one pfam entry
annotations_df = pd.merge(annotations_df, new_df, how='outer')


del enzyme_pfams_list, new_df

    # Kumar 07/08/2022
    # PFAM - RR
    # reaction_id, uniprot_id, Pfams, KO, rhea_id_reaction, kegg_id_reaction, rhea_confirmation, kegg_confirmation,
    #  KO_prediction

print('Loading pfam-RR annotations...')
pfam_rules_df = pd.read_csv("pfam_RR_annotation_file.csv", index_col=None)
pfam_rules_df = pfam_rules_df.rename(columns={pfam_rules_df.columns[1]: 'uniprot_id', pfam_rules_df.columns[2]: 'uniprot_enzyme_pfams_acc'})
pfam_rules_df['reaction_id'] = pfam_rules_df.reaction_id.astype('str')

    # filter type of anotations (strict, medium, loose)
#    if Options.pfam_RR_annotation == 'strict':
#        pfam_rules_df = pfam_rules_df[pfam_rules_df.experimentally_validated]
#    elif Options.pfam_RR_annotation == 'medium':
#        pfam_rules_df = pfam_rules_df[pfam_rules_df.experimentally_validated | pfam_rules_df.Pfam_ec_prediction]
#    else:  # loose
#        pass  # they are all there

        # convert pfam_acc to pfam
pfam_dict = pd.read_csv(pfam_dict_file, index_col=None)
pfam_dict.index = pfam_dict.Acc.apply(lambda x: x.split('.')[0])  # Acc,Name,Desc


    # Kumar
    # Check if the seprataer is working in this case
    # uniprot_enzyme_pfams_acc_list = pfam_rules_df.uniprot_enzyme_pfams_acc.apply(gizmos.pd_to_list, separator=' ')
uniprot_enzyme_pfams_acc_list = pfam_rules_df.uniprot_enzyme_pfams_acc.apply(gizmos.pd_to_list, separator=';')

pfam_rules_df['uniprot_enzyme_pfams_list'] = [[k for k in row if k in pfam_dict.index] for row in uniprot_enzyme_pfams_acc_list]
pfam_rules_df['uniprot_enzyme_pfams'] = pfam_rules_df.uniprot_enzyme_pfams_list.apply(';'.join)

    # Expand df so there is only one pfam per row.
lens = [len(item) for item in pfam_rules_df.uniprot_enzyme_pfams_list]
pfams_flat = [n for sub in pfam_rules_df.uniprot_enzyme_pfams_list for n in sub]

new_df = pd.DataFrame({'uniprot_id': np.repeat(pfam_rules_df.uniprot_id, lens), 'pfam_rule': pfams_flat}).drop_duplicates()

pfam_rules_df = pd.merge(pfam_rules_df, new_df, how='outer')  # on uniprot_id.
del pfams_flat, uniprot_enzyme_pfams_acc_list, new_df
del pfam_rules_df['uniprot_enzyme_pfams_list'], pfam_rules_df['uniprot_enzyme_pfams_acc']


print('Integrating gene annotations with reaction rules...')

merged_df = pd.merge(annotations_df, pfam_rules_df, how='inner')    # on pfam_rule

del merged_df['pfam_rule']
        # now each rule has suspect genes


print('Integrating correlations with gene annotations and reaction rules...')
merged_df = pd.merge(merged_df, correlation_df, how='inner', on='gene')    # on gene

print('Duplicate cleanup...')
merged_df = merged_df.drop_duplicates()         # annotations_df merge produces duplicates due to pfam_rule

print("##### This is from LOAD_ENZYE_INPUT #####")
print(merged_df)

print('Loading transitions...')
transitions_df = pd.read_csv("Jeon_tomato/Fungal/Fungal_pathMassTransition.csv", index_col=None, dtype={'reaction_id': str,
                                                                                          'substrate_id': str,
                                                                                          'product_id': str})
transitions_df['reaction_substrate'] = transitions_df.reaction_id + '_' + transitions_df.substrate_id

transitions_df = transitions_df.rename(columns={'substrate': 'ms_substrate', 'product': 'ms_product',
                                                        'substrate_mnx_id': 'substrate_mnx',
                                                        'substrate_mm': 'substrate_mnx_mm',
                                                        'product_mnx_id': 'product_mnx',
                                                        'product_mm': 'product_mnx_mm'})

transitions_df = transitions_df.rename(columns={'mass_transition_round': 'expected_mass_transition'})
del transitions_df['mass_transition']
    
    # RR db (validated rules)
    # reaction_id, substrate_id, diameter, direction, smarts_id, rxn_smarts, validated
    # ++ rule_substrate_str, rule_product_str, reaction_substrate
print('Loading rules...')
rules_df = pd.read_csv("test_db/format_database/ValidateRulesWithOrigins.csv", index_col=None, dtype={'reaction_id': str, 'substrate_id': str, 'validated': bool})

    # Kumar
    # select only validated reactions (118913/138832)
rules_df = rules_df[rules_df.validated]

rules_df['rule_substrate_str'] = rules_df.rxn_smarts.apply(gizmos.get_subs_string)
rules_df['rule_product_str'] = rules_df.rxn_smarts.apply(gizmos.get_prod_string)
rules_df['reaction_substrate'] = rules_df.reaction_id + '_' + rules_df.substrate_id
rules_df = rules_df.rename(columns={'rule_substrate_str': 'RR_substrate_smarts', 'rule_product_str': 'RR_product_smarts'})

    # MOL
df_to_mol = rules_df[['smarts_id', 'rxn_smarts']].drop_duplicates().copy()
df_to_mol['rxn'] = df_to_mol.rxn_smarts.apply(Chem.ReactionFromSmarts)
rules_df = pd.merge(rules_df, df_to_mol)
del df_to_mol


    # BASE RULES
print('Getting base rules...')
map_df = pd.read_csv("test_db/format_database/base_rules.csv", index_col=None, dtype={'is_base': bool})
    # # # convert to set
map_df['smarts_is_in'] = map_df.smarts_is_in.apply(gizmos.pd_to_set, separator=';')
map_df['smarts_has'] = map_df.smarts_has.apply(gizmos.pd_to_set, separator=';')
map_df['identity'] = map_df.identity.apply(gizmos.pd_to_set, separator=';')
print(map_df.is_base)
base_smarts_id = map_df.smarts_id[map_df.is_base]
print(base_smarts_id)

base_rules_df = rules_df[rules_df.smarts_id.isin(base_smarts_id)][['smarts_id', 'rxn_smarts']].copy()
base_rules_df = base_rules_df.drop_duplicates()
base_rules_df['rxn'] = base_rules_df.rxn_smarts.apply(Chem.ReactionFromSmarts)

    # MERGE
print('Integrating rules and transitions...')
rt_df = pd.merge(rules_df, transitions_df, how='inner')    # on reaction_id, substrate_id


print("##### This is from load and merge transition ####")
print("#### RT_DF ####")
print(rt_df)
print("#### base_rules_df ####")
print(base_rules_df)

merged_df = merged_df.loc[merged_df['correlation'] > 0.70]
print(merged_df)

enzyme_df_subs = merged_df.rename(columns={'ms_name': 'ms_substrate', 'correlation': 'correlation_substrate', 'P': 'P_substrate'})
enzyme_df_prod = merged_df.rename(columns={'ms_name': 'ms_substrate', 'correlation': 'correlation_product', 'P': 'P_product'})

        # merge on reaction_id, gene, and ms_substrate/product
subs_corr_df = pd.merge(rt_df, enzyme_df_subs, on='reaction_id', how='inner')
prod_corr_df = pd.merge(rt_df, enzyme_df_prod, on='reaction_id', how='inner')

        # Here, results may include same rule-gene-coexp data, but through different RR_enzyme annotation
rt_df = pd.merge(subs_corr_df, prod_corr_df, how='outer')  # outer allows for unilateral coexpression
print(rt_df)
        
        # Now we use the allowed transitions to clean the map by removing requirements that are not in the allowed list
        # for this, we remove smarts_id from smarts_has
allowed_smarts_id = set(rt_df.smarts_id)
map_df['smarts_has'] = map_df.smarts_has.apply(lambda x: x.intersection(allowed_smarts_id))

print(map_df)