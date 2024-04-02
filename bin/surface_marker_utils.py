import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import stats
import random
import os
import seaborn as sns
import csv


'''
This file contains a number of minor helper functions.
'''


def flatten(l):
    return [item for sublist in l for item in sublist]

def list_to_csv(listname,filename):
    with open(filename, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
    for val in listname:
        writer.writerow([val])
        
def gaussian(x,mu,b):
    return np.exp(-(x-mu)**2/(2*b**2))/(b*np.sqrt(2*np.pi))

def kde(x,width):
    x = np.array(x)
    xx = np.linspace(min(x),max(x),1000);
    y = np.zeros(1000)
    for i in range(0,x.size):
        point = x[i];
        y += gaussian(xx, point, width)
    riemann = sp.integrate.simps(y,xx)
    y = y/riemann;
    return [xx,y]
    
def plot_tau_kde(val_list, title, path, width = 0.025):
    plt.figure(figsize = (5,5))
    [xx,yy] = kde(val_list, width);
    plt.plot(xx,yy, linewidth= 5);
    plt.title(title)
    plt.ylabel('Probability Density')
    plt.xlabel('Tau Score')
    full_path = '../results/' + path;
    plt.savefig(full_path)
    
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def heatmap_plot(mat, title_str, name, organ_type, brain_cell_type, brain_cutoffs, cell_cutoffs):
    plt.figure(figsize = (10,10));
    sns.heatmap(mat);
    plt.title(title_str);
    plt.xlabel('GTEx {} Tau cutoff'.format(organ_type));
    plt.ylabel('{}-specific Tau cutoff'.format(brain_cell_type));
    plt.xticks(np.arange(0,brain_cutoffs.size), brain_cutoffs);
    plt.yticks(np.arange(0,cell_cutoffs.size), cell_cutoffs);
    full_path = '../results/' + name + '.svg'
    plt.savefig(full_path, transparent =True, dpi = 600)
    

def uniprot_conversion(gene_list, g2u, oldg2u):
    genes_u = []
    for gene in gene_list:
        if gene in g2u:
            genes_u.append(g2u[gene])
        elif gene in oldg2u:
            genes_u.append(oldg2u[gene])
        else:
            genes_u.append('null')
    return genes_u
