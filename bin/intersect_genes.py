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
This function makes a 2d histogram from a given x and y.
'''


def myplot(x, y, s=16, bins=100):
        heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    #     heatmap = gaussian_filter(heatmap, sigma=s)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        return heatmap.T, extent
   


    
'''
This function makes a scatterplot of tau values for genes falling into different categories.
'''

def scatterplot_tau(all_data, heatmap2d = False, to_denote = []):
    
    print('''Creating scatterplots of tau-score with all genes or only genes expressed most in brain/{}
              out of all tissues/cell-types, and with all types of proteins
              or only transmembrane proteins...'''.format(all_data['brain_cell_type']))
    plt.rcParams.update({'font.size': 22})

    from scipy.ndimage.filters import gaussian_filter
    import matplotlib.cm as cm
    
    
    
    g2u = all_data['g2u']
    oldg2u = all_data['oldg2u']
    
    #creates unified dataframe to store both values together and makes scatterplot
    for uniprot_list, name in zip([all_data['full_uniprots'], 
                                   all_data['max_uniprots'], 
                                   all_data['full_uniprots_TMU'], 
                                   all_data['max_uniprots_TMU']],
                                  ['All', 'Maximum', 
                                   'All_TMU', 'Maximum_TMU']):

        #empty matrix to store values
        full_taus = np.zeros((len(uniprot_list), 2));
        plt.figure(figsize = (10,10)) 
        
        #uniprot conversions for genes to denote
        denoted_u = []
        denoted_g = []
        for g in to_denote:
            if g in g2u:
                denoted_u.append(g2u[g])
                denoted_g.append(g)
            elif g in oldg2u:
                denoted_u.append(oldg2u[g])
                denoted_g.append(g)
            else:
                print('Can\'t denote {} in scatterplot, no matching uniprot for this name'.format(g))
        #iterate through all uniprot IDs, grab the tau value from each dataset, and store appropriately
        
        for i, uniprot in tqdm(enumerate(uniprot_list)):
            full_taus[i,0] = all_data['tau_dfs'][0].iloc[all_data['u2idx_dicts'][0][uniprot]]['Tau']
            full_taus[i,1] = all_data['tau_dfs'][1].iloc[all_data['u2idx_dicts'][1][uniprot]]['Tau']
            
            if uniprot in denoted_u:
                gene_name = denoted_g[denoted_u.index(uniprot)]
                plt.scatter(full_taus[i,1],full_taus[i,0], c = 'k', 
                            s = 50, alpha= 1)
                plt.annotate(gene_name, (full_taus[i,1]-0.08,full_taus[i,0]+0.02),
                             c = 'k')


        #creating dataframe, dropping any zeros
        full_tau_df = pd.DataFrame(full_taus, columns = ['GTEx', 'Brain_RNA-Seq' ])
        full_tau_df = full_tau_df.loc[~(full_tau_df==0).any(axis=1)]

        #full scatterplot of all common genes
        qualifier = ''
        gene_type = ''
        beginner = ''
        if 'All' in name:
            size = 75
            beginner = 'All'
        if 'Maximum' in name:
            size = 100
            qualifier = 'with maximum expression'
        if 'TMU' in name:
            gene_type = 'Transmembrane Protein'
            
        if heatmap2d:
            heatmap, extent = myplot(full_tau_df['Brain_RNA-Seq'],
                                     full_tau_df['GTEx'], s= 0.8, bins = 60)
            plt.imshow(heatmap, extent = extent, origin = 'lower', cmap = cm.Blues)
            plt.colorbar(shrink = 0.7)
        else:
            sns.scatterplot(x = 'Brain_RNA-Seq', y = 'GTEx', alpha = 0.3*np.log(231) / np.log(full_tau_df.shape[0]), s= size, 
                            data = full_tau_df)
            
        plt.ylabel('{} Specificity (GTEx Tau Score)'.format(all_data['organ_type']))
        plt.xlabel('{} Specificity (Brain RNA-Seq Tau Score)'.format(all_data['brain_cell_type']))
        
        plt.title('{} {} genes {} in {} and {} ({})'.format(beginner,
                                                            gene_type, 
                                                            qualifier,
                                                            all_data['organ_type'],
                                                         all_data['brain_cell_type'],
                                                            len(uniprot_list)))
        path = '../results/tau_scatterplot_' + name + '.pdf'
        plt.savefig(path)

    print('Done!')
    
    return all_data


'''
This function plots how many genes would be left after different tau cutoffs. Default is minimum 0.6 and maximum 0.9.
'''

def cutoff_heatmap(all_data, tau_min=0.6, tau_max=0.9):
    #plotting number of genes left after different cutoff values between 0.6 and 0.9
    s = 9;
    cutoff_mat = np.zeros((s,s));
    organ_cutoffs = np.linspace(0.6,0.9, s)
    cell_cutoffs = np.linspace(0.9,0.6, s)
    
    for m, organ_cutoff in enumerate(organ_cutoffs):
        for n, cell_cutoff in enumerate(cell_cutoffs):
            cutoff_mat[m, n] = sum([gtex > organ_cutoff and brs > cell_cutoff for 
                                    (gtex, brs) in zip(all_data['master_tau_df']['GTEx'], 
                                    all_data['master_tau_df']['Brain_RNA-Seq'])])

    plt.rcParams.update({'font.size': 10})
    
    print('There are {} transmembrane proteins with maximum in {}/{}.\n'.format(all_data['master_tau_df'].shape[0], all_data['organ_type'], 
 all_data['brain_cell_type']))
    
    heatmap_plot(cutoff_mat[::-1,::-1], 'Number of Candidate Markers with Different Tau Cutoffs','cutoff_heatmap', 
                  all_data['organ_type'], all_data['brain_cell_type'], organ_cutoffs, cell_cutoffs)


'''  
This function takes master tau df and subsets to the genes which pass the tau cutoff. 
It saves those genes into a 'specific_tau' dataframe.
'''

def get_specific_tau(all_data, organ_tau_cutoff = 0.7,
                          cell_tau_cutoff = 0.7, cluster_tau = False):
#     import sys

#     if not sys.warnoptions:
#         import warnings
#         warnings.simplefilter("ignore")

    # creating new dataframe to store only genes decided as specific

   
    # if clustering is used
    if cluster_tau:
        pred_clusters = KMeans(n_clusters = n_clusters).fit_predict(all_data['master_tau_df'].values[:,:2])
        all_data['master_tau_df']['Cluster'] = pred_clusters;

        #plotting clustering results
        plt.figure(figsize = (10,10))
        sns.scatterplot(x = 'Brain_RNA-Seq', y = 'GTEx', hue = 'Cluster',
                        palette = 'Set2', data = all_data['master_tau_df']);
        plt.title('Tau values')
        plt.savefig('tau_clustering.pdf')

        #finding "best cluster" which will we say to have maximum average tau values
        best_cluster = 0
        best_score = -1;
        for cluster in range(0,n_clusters):
            score = sum(np.ravel(all_data['master_tau_df'][all_data['master_tau_df']['Cluster'] == cluster].values[:,:4]))
            if score > best_score:
                best_cluster = cluster;
                best_score = score;

        #creating new dataframe to store just this cluster
        specific_tau = all_data['master_tau_df'].loc[all_data['master_tau_df'].Cluster == best_cluster]

    # otherwise just use individual cutoffs for GTEx and Brain RNA-seq
    else:
        rows = []
        for i in range(0, all_data['master_tau_df'].shape[0]):
            row = all_data['master_tau_df'].iloc[i]
            if row['GTEx'] > organ_tau_cutoff and row['Brain_RNA-Seq'] > cell_tau_cutoff:
                rows.append(i)
        specific_tau = all_data['master_tau_df'].iloc[rows,:]

    all_data['specific_tau'] = specific_tau.copy()
    
    
    return all_data
    

''' 
This function takes the genes in specific tau df and intersects with mass spec data to produce a prioritized candidate list.
It saves this list as a csv and as a dataframe called 'candidates'.
'''

def get_final_list(all_data, exp_name, topk= 20, use_scimilarity = False,
                  show_heatmap = False):
    
    if use_scimilarity:
        specific_tau = all_data['scRNA_data'][all_data['scRNA_data']['attribution'] > 0.001]
        specific_uniprots = []
        for g in specific_tau['gene']:
            if g in all_data['g2u']:
                specific_uniprots.append(all_data['g2u'][g])
            elif g in all_data['oldg2u']:
                specific_uniprots.append(all_data['oldg2u'][g])
            else:
                specific_uniprots.append('null')
        specific_tau= specific_tau.set_index('gene')
    else:
        specific_tau = all_data['specific_tau'].copy()
        specific_uniprots = [u for (u,g) in np.ravel(specific_tau.index)]
    #isolating uniprot IDs for just specific genes
   
    
        # adding in graham's rank and attribution
        g2sc_rank = {g : r+1 for r,g in enumerate(all_data['scRNA_data']['gene'])}
        g2sc_attr = {g : a for g,a in zip(all_data['scRNA_data']['gene'], all_data['scRNA_data']['attribution'])}
        ranks = []; attributions = []
        for u in specific_uniprots:
            if u in all_data['u2g']:
                g = all_data['u2g'][u]
            elif u in all_data['oldu2g']:
                g = all_data['oldu2g'][u]
            else:
                g = 'null'

            if g in g2sc_rank:
                ranks.append(g2sc_rank[all_data['u2g'][u]])
                attributions.append(g2sc_attr[all_data['u2g'][u]])
            else:
                ranks.append(23000)
                attributions.append(0)
        specific_tau['scRNAseq Rank'] = ranks ; specific_tau['scRNAseq Attribution'] = attributions;
    

    #searching mass_spec data
    for name, data in zip(all_data['mass_spec_names'], all_data['mass_spec']):
        specific_tau.loc[:,name] = [int(uniprot in data) for uniprot in specific_uniprots]

    #calculate sum of appearances and sort
    specific_tau.loc[:,'Sum in Mass Spec'] = np.sum(specific_tau.values[:,-len(all_data['mass_spec_names']):],
                                                    axis = 1);
    #add in column for whether it's a transmembrane protein
    specific_tau.loc[:,'Transmembrane'] = [u in all_data['allTMU'] for u in specific_uniprots]
    
    #sort by a specific column
    specific_tau = specific_tau.sort_values('Sum in Mass Spec')
    
    #saves csv
    path = '../results/{}_{}.csv'.format(exp_name, all_data['brain_cell_type'])
    specific_tau.to_csv(path)

    #plot a heatmap of the top 20 candidates by appearance in mass-spec. this should be everything except last 2 columns
    if show_heatmap:
        plt.figure(figsize = (6,6));
        sns.heatmap(specific_tau.iloc[-topk:,-(len(all_data['mass_spec_names'])+2):-2]);
        plt.title('Top 20 Candidate Markers, Ordered by Appearance in Mass-Spec Data');

    all_data['candidates'] = specific_tau.copy()
    
    return all_data


'''
This function takes in just one gene and returns information about that gene's expression across GTEx, HPA, and Brain-RNA-seq,
as well as the appearance of that gene in different mass-spectrometry datasets.
'''

def single_gene_info(g, all_data):
    
    #get the uniprot ID and maybe an alternate gene name
    if g in all_data['g2u']:
        u = all_data['g2u'][g]
    elif g in all_data['oldg2u']:
        u = all_data['oldg2u'][g]
    else:
        print('Could not find a uniprot associated with this gene!')

    
    if u in all_data['u2g']:
        if all_data['u2g'][u] != g:
            alt_g = all_data['u2g'][u]
        else:
            alt_g = 'None'
    elif u in all_data['oldu2g']:
        if all_data['oldu2g'][u] != g:
            alt_g = all_data['oldu2g'][u]
        else:
            alt_g = 'None'
    else:
        alt_g = 'None'

    not_found_gene_exp = []
    
    
    #iterates through the gene expression dataframes to get the expression profiles
    for u2idx_dict, df, tau_df, df_name, gene_list in zip(all_data['u2idx_dicts'],
                                                 all_data['gene_exp_dfs'],
                                                all_data['tau_dfs'],
                                                 all_data['gene_exp_df_names'],
                                              all_data['gene_exp_gene_lists']):
        gene_list = list(gene_list)
        if g in gene_list:
            gene_exp_idx = gene_list.index(g)
        elif alt_g in gene_list:
            gene_exp_idx = gene_list.index(alt_g)
        else:
            not_found_gene_exp.append(df_name)
            continue
        profile= list(df.iloc[gene_exp_idx].values)

        if u in u2idx_dict:
            tau_idx = u2idx_dict[u]
        tau = tau_df.iloc[tau_idx]['Tau']
         
        # making bar graphs for the gene expression
        plt.figure(figsize = ((len(profile)/2, 3)))
        plt.title('{} gene expression in {}: Tau score of {:.3f}'.format(g, df_name, tau))

        plt.bar(x= np.arange(len(profile)), height = profile)
        plt.xticks(np.arange(len(profile)), df.columns, rotation = 90)  
        
    if len(not_found_gene_exp) > 0:
        print('Gene was not found in : {}'.format(not_found_gene_exp))
    
    
    #looking at the mass-spec appearance of the gene and plotting it
    mass_spec_mat = np.zeros((len(all_data['mass_spec_names']), 1))
    for i, (name, data) in enumerate(zip(all_data['mass_spec_names'], all_data['mass_spec'])):
        mass_spec_mat[i] = int(u in data) 
    plt.figure(figsize = (0.4,3))
    sns.heatmap(mass_spec_mat)
    plt.yticks(np.arange(len(all_data['mass_spec_names']))+0.5, all_data['mass_spec_names'], rotation = 0);




