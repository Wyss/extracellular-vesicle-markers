#importing libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy as sp
from scipy import stats
import os
import seaborn as sns
import csv
from tqdm import tqdm
from sklearn.cluster import KMeans
from itertools import compress

#custom functions to import proteomics/expression/annotation data
from import_surfacemarker_data import import_proteins, import_surfaceome, import_GTEx, import_brain_RNAseq, import_hpa_data, import_MS_DT, import_IPSC_data, import_HPPP, import_past_CSF_MS


'''

This function sets up all the datasets and compiles them into the master 'all_data' dictionary to make them easily accessible.

'''
def setup_datasets(brain_cell_type = 'neurons', organ_type= 'Brain'):
    print('Importing protein annotations/ID conversions...')

    oldTMU, TMU, ec_tmu, i2u, g2u, u2g, oldg2u, oldu2g = import_proteins();
    allTMU = set(oldTMU + TMU)

    surfaceome_TMU = import_surfaceome();

    print('Done!')
    
    #GTEx data
    compiled_gtex, gtex_genes = import_GTEx(); 

    #Brain RNA seq data
    compiled_brs, brs_genes = import_brain_RNAseq() 

    #HPA data
    hpa_ihc, hpa_ihc_genes, hpa_rna, hpa_rna_genes, hpa_consensus, hpa_consensus_genes  = import_hpa_data()
    
    scrna_data = pd.read_csv('../data/scimilarity_attributions/{}_gene_attributions.csv'.format(brain_cell_type))
    
    
    print('Importing mass spectrometry data...')

    #captocore-isolated exosome proteomics 
    mass_spec_plasma_cc, mass_spec_csf_cc, mass_spec_neuron_cc = import_MS_DT();

    #iPSC proteomics
    IPSC_uniprots = import_IPSC_data(brain_cell_type, g2u, oldg2u)

    #HPPP data
    HPPP_plasma_cc, HPPP_ev_cc = import_HPPP();

    #Past CSF total/EV datasets
    CSF_totals, CSF_EVs, CSF_total_names, CSF_EV_names = import_past_CSF_MS(i2u)
    mass_spec_names = ['DT Mass Spec - Plasma', 'DT Mass Spec - CSF',
                       'DT Mass Spec - Neuron Culture',
                       'HPPP - Plasma', 'HPPP - CSF', 'iPSC Proteomics'] + CSF_total_names + CSF_EV_names;

    mass_spec_data = [mass_spec_plasma_cc, mass_spec_csf_cc, mass_spec_neuron_cc, 
                     HPPP_plasma_cc, HPPP_ev_cc, IPSC_uniprots] + CSF_totals + CSF_EVs
    print('Done!')
    
    
    
    all_data= {'g2u': g2u,
               'oldg2u': oldg2u,
               'u2g': u2g,
               'oldu2g': oldu2g,
               'TMU': ec_tmu,
               'allTMU': allTMU,
               'surfaceome': surfaceome_TMU,
               'GTEX_data': compiled_gtex,
               'GTEX_names': gtex_genes,
               'BRS_data': compiled_brs,
               'BRS_names': brs_genes,
               'HPA_IHC_data': hpa_ihc,
               'HPA_IHC_names': hpa_ihc_genes,
               'HPA_RNA_data': hpa_rna,
               'HPA_RNA_names': hpa_rna_genes,
               'HPA_consensus_data': hpa_consensus,
               'HPA_consensus_names': hpa_consensus_genes,
               'gene_exp_dfs': [compiled_gtex, compiled_brs, hpa_rna, hpa_ihc, 
                                hpa_consensus],
               'gene_exp_gene_lists': [gtex_genes, brs_genes, hpa_rna_genes, 
                                       hpa_ihc_genes, hpa_consensus_genes],
               'gene_exp_df_names' : ['GTEx', 'Brain_RNA-Seq', 'HPA_RNA', 'HPA_IHC', 
                                      'HPA_Consensus'],
               'scRNA_data': scrna_data,
               'mass_spec': mass_spec_data,
               'mass_spec_names': mass_spec_names,
               'iPSC_uniprots': IPSC_uniprots,
               'brain_cell_type': brain_cell_type,
               'organ_type': organ_type
              }
    return all_data

               
               
    
    
    
