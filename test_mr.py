#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from multiprocessing import Pool, Manager
import sqlite3
from multiprocessing.shared_memory import SharedMemory
from multiprocessing.managers import SharedMemoryManager

# project specific module
import mutual_ranks
import gizmos


######
conn = sqlite3.connect("Jeon_bac.sqlite")
cur = conn.cursor()
sql1 = "select * from jeon_bac_diff_exp_MR_weights_DR_10 where metabolite in ('M723T804') and gene in ('Solyc01g108600.2', 'Solyc02g093680.2', 'Solyc07g005760.2');"


#Solyc07g062520  Cluster 899 DR 5
#M544T534 Solyc07g008570.2 Solyc11g013110.1 Solyc07g062520.2 Solyc11g010960.1
#Solyc07g043460  Cluster 1522 DR 10
#M266T597 Solyc03g043740.2 Solyc03g121180.2 Solyc07g043460.2
#Solyc12g006460 Cluster 2114 DR 10
#M360T818 Solyc05g009080.2 Solyc06g019200.2 Solyc12g006460.1 Solyc03g117690.2
#Solyc07g043420 Cluster 468 DR 10
#M266T597 Solyc07g043420.2 Solyc03g043740.2 Solyc03g121180.2 Solyc02g086180.2
#Solyc12g006470  Cluster 276 DR 10
#M462T635 M492T803 M330T805 M484T775 Solyc12g100120.1 Solyc04g070980.2 Solyc06g054320.1 Solyc07g045540.2 Solyc12g006470.1 Solyc01g106360.2 Solyc04g056390.2 Solyc07g043480.1 Solyc03g025350.2 Solyc08g015940.2 Solyc06g083690.2
#Solyc07g043490  Cluster 357 DR 10
#M266T805 Solyc02g036470.2 M462T635 M492T803 Solyc01g068030.2 Solyc12g096510.1 Solyc08g061100.2 M484T775 Solyc04g009870.1 Solyc02g082920.2 Solyc04g070980.2 Solyc06g054320.1 Solyc01g105360.2 Solyc08g080210.2 Solyc07g063530.2 Solyc05g053100.2 Solyc05g018650.2 Solyc09g082730.2 Solyc04g080390.2 Solyc08g078210.2 Solyc07g043490.1 Solyc12g009440.1 Solyc08g076840.2 Solyc02g068050.2 Solyc03g071720.2 Solyc09g098040.2 Solyc04g064720.2 Solyc08g077060.2 Solyc06g005230.1 Solyc03g025350.2 Solyc01g008790.2 Solyc06g075680.1 Solyc09g065920.2 Solyc02g083280.2 Solyc03g121270.2 Solyc03g005340.2 Solyc02g093320.2 Solyc10g076250.1
#Solyc07g043480 Cluster 576 DR 10
#M462T635 M492T803 M484T775 Solyc04g070980.2 Solyc06g054320.1 Solyc01g098560.2 Solyc07g043480.1 Solyc03g025350.2
#Solyc07g043500 Cluster 375 DR 10
#M382T711 M342T713 Solyc07g043500.1 Solyc07g007250.2 Solyc04g081890.1 Solyc07g032490.2 Solyc09g098590.2
#GAME 9 (regulator of SGA Solyc01g090340) (new) cluster 1606 DR 5; Cluster 1140 DR 10; Cluster 728 DR 25; Cluster 665 DR 50; Cluster 673 DR 100
#M190T392 Solyc01g096370.2 Solyc08g005130.1 Solyc01g090340.2 Solyc02g062950.2
#
#Chlorogenic acid biosynthetic pathway 
#Bacterial (new )Solyc07g005760  Cluster 3165 DR 10
#M723T804 Solyc01g108600.2 Solyc02g093680.2 Solyc07g005760.2


query = cur.execute(sql1)
results_df = pd.DataFrame.from_records(data = query.fetchall(), columns=["metabolite", "gene", "MR", "weight"])
print(results_df)
results_df.to_csv("HQT_Cluster_3165_DR_10.csv", index=False)
conn.close()
#######

#metabolites = [i for i in list(results_df['metabolite'].squeeze()) if i is not None]
#genes = [i for i in list(results_df['gene'].squeeze()) if i is not None]
#all_objects = metabolites + genes
#all_objects = list(dict.fromkeys(all_objects))
#
#edges_df = mutual_ranks.get_MR_from_correlations(results_df, all_objects)
#print(edges_df)
#        
# Kumar
# Decay rate could be any value but generally it is either {5, 10, 25, 50, 100} from the Wisecaver (2017) paper.    
#weights_df = mutual_ranks.get_weight_from_MR(edges_df, decay_rate = 25)
#print(weights_df)