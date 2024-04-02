from setup import setup_datasets
from tau_calculations import make_taus
from tau_calculations import count_gene_categories
from tau_calculations import make_unified_gene_list,unify_tau_dataframes
from intersect_genes import scatterplot_tau, cutoff_heatmap
from intersect_genes import get_specific_tau, get_final_list
from scipy.ndimage import gaussian_filter
import matplotlib.cm as cm
from intersect_genes import myplot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np




'''
S3 will unify based on genes which are maximally expressed in both the brain and the cell-type of choice.
Then it will use the tau-cutoffs of 0.7 for each to get specific genes and write a csv with those gene names.
'''

def s3(all_data, brain_cell_type):
    all_data = unify_tau_dataframes(all_data, all_data['max_uniprots'])
    all_data = get_specific_tau(all_data, 
                                     organ_tau_cutoff = 0.7,
                                      cell_tau_cutoff = 0.7,
                                     cluster_tau = False)

    f = open('../supp_tables/table_s3_specific_brain_{}.csv'.format(brain_cell_type), 'w+')
    f.write('Number,Gene Name, Transmembrane Domain?\n')

    for i, (u,g) in zip(np.arange(all_data['specific_tau'].shape[0]), all_data['specific_tau'].index):
        if u in all_data['u2g']:
            f.write('{},{},{}\n'.format(i+1,all_data['u2g'][u], u in all_data['allTMU']))
        elif u in all_data['oldu2g2g']:
            f.write('{},{},{}\n'.format(i+1,all_data['oldu2g'][u], u in all_data['allTMU']))
    f.close()
    
'''
S4 will find the sCimilarity scores for all the genes in each brain cell-type and write a csv which has them
and whether they have a transmembrane domain or not.
'''
    
def s4(all_data, brain_cell_type):
    f = open('../supp_tables/table_s4_scimilarity_{}.csv'.format(brain_cell_type), 'w+')
    f.write('Number,Gene Name, SCimilarity Attribution Score, Transmembrane Domain?\n')
    TMU_class = []
    for i, g, score in zip(np.arange(all_data['scRNA_data'].shape[0]), 
                           all_data['scRNA_data']['gene'],
                           all_data['scRNA_data']['attribution']):
        if g in all_data['oldg2u']:
            tm = all_data['oldg2u'][g] in all_data['allTMU']
            f.write('{},{},{},{}\n'.format(i+1,g, score, tm))
        elif g in all_data['g2u']:
            tm = all_data['g2u'][g] in all_data['allTMU']
            f.write('{},{},{},{}\n'.format(i+1,g, score, tm))
    f.close()
    
'''
S5 will unify based on genes which are maximally expressed in both the brain and the cell-type of choice.
Then it will subset to the specific genes using a 0.7 tau score cutoff. Next, it will subset this list to only those
with scRNAseq attribution > 0.001 for the given cell-type and write a csv with those genes.
'''

def s5(all_data, brain_cell_type):
    all_data = unify_tau_dataframes(all_data, all_data['max_uniprots'])
    all_data = get_specific_tau(all_data, 
                                     organ_tau_cutoff = 0.7,
                                      cell_tau_cutoff = 0.7,
                                     cluster_tau = False)

    all_data = get_final_list(all_data, exp_name = '20230730', topk= 20, use_scimilarity = False)

    table = all_data['candidates'][all_data['candidates']['scRNAseq Attribution'] > 0.001]
    f = open('../supp_tables/table_s5_scimilarity_and_bulk_specific_{}.csv'.format(brain_cell_type), 'w+')
    f.write('Number,Gene Name, Transmembrane Domain?\n')
    uniprots = [u for u,g in table.index]
    for i, u in zip(np.arange(table.shape[0]), uniprots):
            f.write('{},{},{}\n'.format(i+1,all_data['u2g'][u], u in all_data['allTMU']))
    f.close()

'''
S6 will unify based on genes which are maximally expressed in both the brain and the cell-type of choice
and have a transmembrane domain. Then it will subset to the specific genes using a 0.7 tau score cutoff.
Finally, it will save this list of genes along with information about mass-spec appearances.
'''

def s6(all_data, brain_cell_type):
    all_data = unify_tau_dataframes(all_data, all_data['max_uniprots_TMU'])
    all_data = get_specific_tau(all_data, 
                                     organ_tau_cutoff = 0.7,
                                      cell_tau_cutoff = 0.7,
                                     cluster_tau = False)

    all_data = get_final_list(all_data, exp_name = '20230730', topk= 20, use_scimilarity = False)
    all_data['candidates'].to_csv('../supp_tables/table_s6_marker_candidates_{}.csv'.format(brain_cell_type))