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

from import_surfacemarker_data import *
from surface_marker_utils import *


'''
This function calculates tau values for a given expression matrix of genes x organs or genes x cell-types.

For each gene, it saves the Tau value, whether that gene was expressed more in the given 'idx' organ or cell-type than any other,
and whether the gene has a valid uniprot ID. This will matter for downstream processing.

Then it returns a dataframe with this information for all the genes (tau_df),
and a subset dataframe which has only the genes which were most in the desired organ or cell type (max_tau_df).
'''
def tau_calc(exp_dataset, gene_list_uniprot, idx):
    
    tau_df = pd.DataFrame(gene_list_uniprot, columns = ['Gene']);
    taus = []
    max_types = []
    valid_uniprot = []
    valid_max = []

    for g, gene in enumerate(gene_list_uniprot): #iterate through every gene as a uniprot ID
        
        #extract the pseudologged expression profile of the gene across cell-type/tissues
        profile = exp_dataset[g,:] 
        
        T = 0;
        max_type = False;
        
        #if the gene has nonzero expression
        if max(profile) > 0.3: 
            #unit-scale profile
            profile = profile/max(profile); 
            
            #calculate Tau
            T = sum((1-profile)/(profile.size-1));
            
            #save whether maximal expression is in right cell-type/tissue
            max_type = np.argmax(profile) == idx;
            
        taus.append(T);
        max_types.append(max_type);
        
        #checking whether gene is a valid uniprot ID
        valid_uniprot.append(gene != 'null');
        valid_max.append(gene != 'null' and max_type)
        if gene == 'Q5HZK2':
            print(profile)
            
    tau_df['Tau'] = taus;
    tau_df['Maximum'] = max_types;
    tau_df['Valid_Uniprot'] = valid_uniprot;
    tau_df['Valid_Max'] = valid_max
    
    #creating separate dataframe for only valid genes with max in the right cell-type/tissue
    max_taus = tau_df.loc[tau_df.Valid_Max,:]
    return tau_df, max_taus



'''
This function processes a given expression dataset.
First it pseudologs the data (pl_data), then figures out which index matches the desired organ or cell-type (type_idx).
Then it uses this data with the tau_calc function to make tau dataframes.
It optionally will show a kernel-density estimate of the tau value distribution.
'''
def data_tau_process(all_data, expression_data, name, genes, brain_cell_type, organ_type, show_plots = True):
    
    pl_data = np.log2(expression_data.values.astype(float)+1);
    genes_uniprot = list_conversion(genes, all_data['g2u'], all_data['oldg2u']);
    
    if name == 'Brain_RNA-Seq':
        type_idx = list(expression_data.columns).index(brain_cell_type);
    elif organ_type in expression_data.columns:
        type_idx = list(expression_data.columns).index(organ_type);
    elif organ_type.lower() in expression_data.columns:
        type_idx = list(expression_data.columns).index(organ_type.lower());
    else:
        type_idx = 0
        print('Organ not found, going with default')
        
    tau, type_tau = tau_calc(pl_data, genes_uniprot, type_idx)
    
    title= name +' {}-Maximum Tau Distribution'.format(organ_type);
    
    if name == 'Brain_RNA-Seq':
        title= name +' {}-Maximum Tau Distribution'.format(brain_cell_type);
        
    path = name+'_max.pdf'
    if show_plots:
        plot_tau_kde(type_tau.Tau.values,title,path)
    return tau, type_tau

'''
This function goes through every gene expression dataset and does the tau calculation on all of them.
It saves both the tau dataframes and the subset dataframes which have only the maximally-expressed genes in 'all_data'.
'''
def make_taus(all_data, show_plots = True):
    #creating lists of expression datasets, names, and gene lists
    brain_cell_type = all_data['brain_cell_type']
    organ_type = all_data['organ_type']
    tau_dfs = []
    max_tau_dfs = []

    #iterating through and calculating tau
    for df, df_name, gene_list in zip(all_data['gene_exp_dfs'],
                                      all_data['gene_exp_df_names'],
                                      all_data['gene_exp_gene_lists']):
        tau_df, max_tau_df = data_tau_process(all_data, df, df_name, gene_list, brain_cell_type, organ_type, show_plots= show_plots);
        tau_dfs.append(tau_df);
        max_tau_dfs.append(max_tau_df)
    all_data['tau_dfs'] = tau_dfs
    all_data['max_tau_dfs'] = max_tau_dfs
    
    return all_data


'''
This function can be used to get a sense of gene counts across a number of categories.
'''
def count_gene_categories(all_data, brain_cell_type = 'neurons'):
    #plotting dataset counts
    fig, ax = plt.subplots(figsize = (12,5))
    ind = np.arange(0,6)*7;
    
    ax.set_xticks(ind+2, ['Total',
                          'Max in {}/Brain'.format(brain_cell_type),
                          'Map to Uniprot ID', 'Max and Map',
                          'Max and TMU', 'Max and Surfaceome TMU'],
                          rotation= 45 )
    
    
    #matrix to store counts
    count_mat = np.zeros((6, len(all_data['tau_dfs'])))
    
    
    for i, (df, max_tau_df, name) in enumerate(zip(all_data['tau_dfs'],
                                                  all_data['max_tau_dfs'],
                                                  all_data['gene_exp_df_names'])): 
        
        count_mat[:4,i] = df.shape[0], sum(df.Maximum), sum(df.Valid_Uniprot), sum(df.Valid_Max)
        count_mat[4,i] = len(intersection(max_tau_df.Gene, all_data['TMU']));
        count_mat[5,i] = len(intersection(max_tau_df.Gene, all_data['surfaceome']))

        height = count_mat[:,i];
        x = ind + i;
        ax = plt.bar(x = x, height = height)

        s = [str(int(c)) for c in count_mat[:,i]]
        y = height*1.2
        for xi, yi, si in zip(x+0.3,y,s):
            plt.text(xi, yi, si, ha= 'right', rotation = 90)

    plt.title('Number of Genes in Different Categories')
    plt.yscale('log');
    plt.ylim([1,10**5.5])
    plt.legend(all_data['gene_exp_df_names']);
    plt.savefig('../results/counts.pdf')
    
'''
This function makes some useful lists of genes that can be used for unification.
Importantly, it also creates a u2idx dictionary that can be used to quickly look up tau-values for a given gene.
'''

def make_unified_gene_list(all_data):
    print('Creating unified gene lists...')

    #dictionaries for each dataset that store what index every gene corresponds to 
    all_data['u2idx_dicts'] = [dict(zip(df.Gene,np.arange(len(df.Gene)))) for df in all_data['tau_dfs']]

    #lists of genes in common between GTEx and brain RNA seq
    all_data['full_uniprots'] = intersection(list(all_data['tau_dfs'][0].Gene),
                                             list(all_data['tau_dfs'][1].Gene))

    #those which are maximum in both for brain/cell-type
    # max_uniprots = list(gtex_brain_tau.Gene)
    all_data['max_uniprots'] = intersection(list(all_data['max_tau_dfs'][0].Gene),
                                             list(all_data['max_tau_dfs'][1].Gene))

    all_data['full_uniprots_TMU'] = intersection(all_data['full_uniprots'],
                                                 all_data['allTMU'])

    all_data['max_uniprots_TMU'] = intersection(all_data['max_uniprots'],
                                                all_data['allTMU'])
    print('Done!')
    
    return all_data
    
'''
This function uses a list of genes to make a unified tau dataframe. 
Each row will be one gene in the list and each column would be tau values for a given expression dataset.
Unification makes future subsetting of genes much easier.
'''

def unify_tau_dataframes(all_data, master_gene_list):

    tau_mat = np.zeros((len(master_gene_list),5))

    #creating new matrix to store Tau values across expression datasets for TMUs with max in brain/cell-type
    print('Unifying tau values into a single matrix...')
    gene_names= []
    
    #going through each uniprot i, finding index in df j, and setting the tau value at i,j
    for i, uniprot in tqdm(enumerate(master_gene_list)):
        
        if uniprot in all_data['u2g']:
            gene_names.append(all_data['u2g'][uniprot]);
        elif uniprot in all_data['oldu2g']:
            gene_names.append(all_data['oldu2g'][uniprot])
        else:
            gene_names.append('null')
            
        for j, tau_df in enumerate(all_data['tau_dfs']):
            if uniprot in all_data['u2idx_dicts'][j]:        
                tau_mat[i,j] = tau_df.iloc[all_data['u2idx_dicts'][j][uniprot]]['Tau'];
    print('Done!')
    assert(len(gene_names) == len(master_gene_list))

    #making dataframe with uniprot IDs and gene names
    master_tau_df = pd.DataFrame(tau_mat, columns = all_data['gene_exp_df_names'])

    master_tau_df['Uniprot ID'] = master_gene_list;
    master_tau_df['Gene name'] = gene_names;
    master_tau_df = master_tau_df.set_index('Uniprot ID')
    master_tau_df = master_tau_df.set_index('Gene name', append= True)

    #summing tau values
    master_tau_df['Tau Sum'] = np.sum(tau_mat, axis = 1);
    master_tau_df = master_tau_df.sort_values('Tau Sum')
    
    all_data['tau_mat'] = tau_mat
    all_data['master_tau_df'] = master_tau_df
    
    return all_data
