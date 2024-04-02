# gene ID conversions

# load in list of all reviewed human proteins
import pandas as pd
import numpy as np
from tqdm import tqdm


'''
This function imports information about proteins from uniprot. 
We use two different annotation files to account for changes in uniprot. 
It also creates gene name/IPI to uniprot ID dictionaries.

The primary way this codebase handles genes is as uniprot IDs, so making conversions is important.
'''
def import_proteins():
    annotations = pd.read_csv('../data/prot_ann/uniprot-organism_homo_sapiens.tab',  sep='\t')
    human_annotations = annotations[annotations.Organism == 'Homo sapiens (Human)']

    TMU = []
    g2u = {}
    for i in tqdm(range(len(human_annotations))):
        row = human_annotations.iloc[i]
        if row.Status == 'reviewed':
            if str(row['Gene names']) != 'nan':
                g2u[row['Gene names'].split(' ')[0]] = row['Entry'];
            if str(row['Transmembrane']) != 'nan':
                TMU.append(row['Entry']);
                
    annotations = pd.read_csv('../data/prot_ann/uniprotkb_homo_sapiens_AND_model_organi_2023_12_06.tsv',  sep='\t')
    ec_tmu = []
    for i in range(annotations.shape[0]):
        row = annotations.iloc[i]
        domain_ann = str(row['Topological domain']).lower()
        if str(row['Transmembrane']) != 'nan' and 'extracellular' in domain_ann and 'secreted' not in domain_ann:
            ec_tmu.append(row['Entry'])

    u2g = {v: k for k, v in g2u.items()};

    i2u_df= pd.read_csv('../data/prot_ann/i2u.csv', header = None).rename(columns = {0:'IPI', 1:'Uniprot'})
    i2u = {}
    for row in range(0,i2u_df.shape[0]):
        i2u[i2u_df.iloc[row]['IPI']] = str(i2u_df.iloc[row]['Uniprot'])[:6];

    old_allProt = pd.read_csv('../data/prot_ann/uniprot-organism_homo_sapiens_old.tab', sep='\t')
    old_allProt = old_allProt[old_allProt.Status == 'reviewed']
    old_allProt = old_allProt.rename(columns={'Gene names':'GeneNames'})

    # split gene names into list
    old_allProt['GeneNames_list'] = [str(x).replace('/',' ').split(' ') for x in old_allProt.GeneNames]

    # make dictionary to map Uniprot IDs to gene names
    oldg2u = {}
    for i in old_allProt.index:
        gl = old_allProt.loc[i,'GeneNames_list']
        for g in gl:
            oldg2u[g] = old_allProt.loc[i,'Entry']
    oldu2g = {v: k for k, v in oldg2u.items()};
    # filter all human proteins to transmembrane proteins only
    TMProt = old_allProt[pd.notnull(old_allProt['Transmembrane'])]

    # filter out mitochondrial proteins
    TMProt = TMProt.rename(columns={'Subcellular location [CC]':'SubcellularLocation'})
    TMProt = TMProt[np.logical_or(TMProt.SubcellularLocation.str.contains("Mitochon") == False, 
                                  pd.isnull(TMProt.SubcellularLocation))]

    # get list of Uniprot IDs for these proteins
    oldTMU= list(TMProt.Entry)
    
    return oldTMU, TMU, ec_tmu, i2u, g2u, u2g, oldg2u, oldu2g
        
'''
This function imports data from our own mass-spectrometry datasets of plasma, CSF, and IPS-neuron derived EVs.
If any isoform of a given protein is found, we simply denote that this protein was found and drop the isoform annotation.
This is because the RNA-sequencing datasets we use are not isoform resolved.
'''

def import_MS_DT():

    mass_spec_plasma = pd.read_csv('../data/DT_MS/201210X_SAM7799_Plasma_CC_125_mPOP-TK_ALLPlasma_LFQ.csv'); 
    mass_spec_plasma_cc = [u.split('-')[0] for u in mass_spec_plasma['Accession'].values];
    # mass_spec_plasma_cc = [u for u in mass_spec_plasma['Accession'].values];

    mass_spec_csf = pd.read_csv('../data/DT_MS/Dima_CSFsamples_SAM7157-8875_ALL_LFQ.csv');
    mass_spec_csf_cc = [u.split('-')[0] for u in mass_spec_csf['Accession'].values];
    # mass_spec_csf_cc = [u for u in mass_spec_csf['Accession'].values];

    neuron_mass_spec= pd.read_csv('../data/DT_MS/raw_iPSneuron_EV_MS.csv')
    mass_spec_neuron_cc= []
    for i in range(0, neuron_mass_spec.shape[0]):
        if str(neuron_mass_spec.iloc[i]['More than 2 peptides/protein in Neuron only']) != 'nan':
            mass_spec_neuron_cc.append(neuron_mass_spec.iloc[i]['accession_number'].split('-')[0])
        else:
            break;
        
    return mass_spec_plasma_cc, mass_spec_csf_cc, mass_spec_neuron_cc

'''
This function imports HPPP proteomic datasets for plasma and plasma derived EVs.
'''

def import_HPPP():
    HPPP_plasma = pd.read_csv('../data/HPPP/human_plasma_proteome_2021-07.csv');
    HPPP_plasma_cc = HPPP_plasma['biosequence_name'].values;
    HPPP_ev = pd.read_csv('../data/HPPP/human_plasma_ev_proteome_2021-06.csv');
    HPPP_ev_cc = HPPP_ev['biosequence_name'].values;

    return HPPP_plasma_cc, HPPP_ev_cc

'''
This function imports several previously published proteomic datasets for CSF.
'''

def import_past_CSF_MS(i2u):
    CSF_total_paths = ['../data/CSF_Data/CSF_total/Gulbrandsen2014/Gulbrandsen2014_S6.csv',
     '../data/CSF_Data/CSF_total/Zhang2007/Zhang2007_S1.csv',
     '../data/CSF_Data/CSF_total/Macron2018/Macron2018_S1.csv',
    '../data/CSF_Data/CSF_total/Schutzer2010/Schutzer2010_S1.csv', '../data/CSF_Data/CSF_total/Begcevic2016/Begcevic2016uniprots.csv']


    CSF_total_names = ['Guldbrandsen2014', 'Zhang2007', 'Macron2018', 'Schutzer2010', 'Begcevic2016']

    CSF_total_acc = ['Guldbrandsen','IPI\nAccession #', 'Protein Accession Number', 'IPI', 'Accession']

    total_IPI = [False,True,False,True, False];


    CSF_EV_paths = ['../data/CSF_Data/CSF_EV/Thompson2020/Thompson2020_S2.csv','../data/CSF_Data/CSF_EV/Muraoka2020/Muraoka2020_S1.csv',
    '../data/CSF_Data/CSF_EV/Thompson2018/Thompson2018_S1.csv', '../data/CSF_Data/CSF_EV/Chiasserini2014/Chiasserini2014_S1.csv' ];

    CSF_EV_names = ['Thompson2020', 'Muraoka2020', 'Thompson2018', 'Chiasserini2014'];
    CSF_EV_acc = ['Accession', 'Uniprot ID', 'Accession', 'uniprot accession']


    CSF_totals = [];

    CSF_EVs = [];

    for i,path in enumerate(CSF_total_paths):
        df = pd.read_csv(path);

        if total_IPI[i]:
            uniprots = [];
            for IPI in df[CSF_total_acc[i]]:
                if IPI in i2u.keys():
                    if i2u[IPI] != None: 
                        uniprots.append(i2u[IPI].split('.')[0]);
        else:
            uniprots = [uniprot[:6] for uniprot in df[CSF_total_acc[i]]];

        CSF_totals.append(set(uniprots))

    for i,path in enumerate(CSF_EV_paths):
        df = pd.read_csv(path);
        uniprots = [uniprot[:6] for uniprot in df[CSF_EV_acc[i]]];
        if CSF_EV_names[i] == 'Chiasserini2014':
            uniprots = [];
            for j in range(0,df.shape[0]):
                row = df.iloc[j]
                if row['Intensity EV1'] ==0 and row['Intensity EV2'] ==0:
                    continue;
                else:
                    uniprots.append(row[CSF_EV_acc[i]])
        CSF_EVs.append(set(uniprots))

    return CSF_totals, CSF_EVs, CSF_total_names, CSF_EV_names



'''
This function imports all surfaceome classifications. 
'''
def import_surfaceome():
    surfaceome = pd.read_csv('../data/surfaceome_classification.csv')
    surfaceome_TMU = [];
    for row in range(0,surfaceome.shape[0]):
        if surfaceome.iloc[row]['Surfaceome Label'] == 'surface':
            surfaceome_TMU.append(surfaceome.iloc[row]['UniProt accession'])
    return surfaceome_TMU
    
'''
This function imports gene expression data from GTEx.
We average expression across all the regions or tissues in a given organ.
We also exclude expression in the nerve, pituitary, and testis.
'''

def import_GTEx():
    gtex = pd.read_csv('../data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct', sep='\t',header=2)
    compiled_gtex_tissues = pd.DataFrame()
    tissue_counts = {}
    
    #averages different regions in the same tissue
    for t1, tissue1 in enumerate(gtex.columns[2:]):
            region1 = tissue1.split(' - ')[0];
            if region1 not in compiled_gtex_tissues.columns:
                compiled_gtex_tissues[region1] = gtex[tissue1].values;
                tissue_counts[region1] = 1;
            else:
                compiled_gtex_tissues[region1] = (compiled_gtex_tissues[region1]*tissue_counts[region1] + gtex[tissue1].values)/ (tissue_counts[region1]+1);
                tissue_counts[region1] += 1;

    # compiled_GTEX_tissues['Brain'] = np.max(GTEx.values[:,['Brain' in col for col in GTEx.columns]], axis = 1)
    
    drop_tissues = ['Nerve', 'Pituitary', 'Testis']
    compiled_gtex_tissues_adjust = compiled_gtex_tissues.drop(drop_tissues, axis = 1)
    
    return compiled_gtex_tissues_adjust, gtex.Description

    

#--------------------------importing Brain RNA seq data and compiling cell types together


'''
This function imports gene expression data from Brain RNA seq.
We average the expression across all the mature cell-types to produce a single figure for a given cell-type.
'''

def import_brain_RNAseq():
    
    # load in FPKM data from brs et al 2016
    brs_2016 = pd.read_csv('../data/Barres2016TableS4.csv', low_memory=False)
    brs_2016 = brs_2016.drop([0,1])

    # label sample columns
    brs_2016.columns = ['Gene','GBM-PT-astrocytes-1','GBM-PT-astrocytes-2','GBM-PT-astrocytes-3','GBM-PT-astrocytes-4', 'SH-astrocytes-1',
                       'SH-astrocytes-2','SH-astrocytes-3','SH-astrocytes-4','fetal-astrocytes-1','fetal-astrocytes-2','fetal-astrocytes-3','fetal-astrocytes-4',
                       'fetal-astrocytes-5','fetal-astrocytes-6',
                       'astrocytes-1','astrocytes-2','astrocytes-3','astrocytes-4','astrocytes-5','astrocytes-6','astrocytes-7','astrocytes-8',
                       'astrocytes-9', 'astrocytes-10','astrocytes-11', 'astrocytes-12',
                       'neurons',
                       'oligodendrocytes-1','oligodendrocytes-2','oligodendrocytes-3','oligodendrocytes-4','oligodendrocytes-5',
                       'microglia-1','microglia-2','microglia-3',
                       'endothelial-1','endothelial-2',
                       'cortex-1','cortex-2','cortex-3','cortex-4']

    # set data type to numeric
    for sample in brs_2016.columns[1:]:
        brs_2016[sample] = pd.to_numeric(brs_2016[sample])

    brs_2016 = brs_2016.drop([23225, 23226]);

    sep = '_'
    celltype_counts = {}
    compiled_brs_2016 = pd.DataFrame()
    
    #averaging different replicates of cell-types togethe
    for col in brs_2016.columns[1:]:
        if '-' in col:
            newcol= sep.join(col.split('-')[:-1])
        else:
            newcol = col;
        if newcol not in compiled_brs_2016.columns:
            compiled_brs_2016[newcol] = brs_2016[col];
            celltype_counts[newcol] = 1;
        else:
            compiled_brs_2016[newcol] = (compiled_brs_2016[newcol]*celltype_counts[newcol] +brs_2016[col]) /  (celltype_counts[newcol]+1);
            celltype_counts[newcol] += 1;
            
    compiled_brs_2016_adjust = compiled_brs_2016.drop(['fetal_astrocytes', 'GBM_PT_astrocytes',
                                                   'SH_astrocytes', 'cortex'], axis = 1)
    return compiled_brs_2016_adjust, brs_2016.Gene,


'''
This function adjusts gene expression data columns to average together all the different brain regions/tissues.
It also drops the testis.
'''

def adjust_hpa_data(hpa_data, brain_parts = []):
    brain_cols = [] 
    if len(brain_parts) > 0:
        for i,tissue in enumerate(hpa_data.columns):
            if any([brain_part in tissue for brain_part in brain_parts]):
                brain_cols.append(i);
    else:
        for i,tissue in enumerate(hpa_data.columns):
            if 'Brain' in tissue:
                brain_cols.append(i);
                
    average_brain_cols = np.max(hpa_data.values[:,brain_cols], axis = 1); 
    hpa_adjust = hpa_data.drop(hpa_data.columns[brain_cols], axis = 1);          
    hpa_adjust = hpa_adjust.drop(['testis'], axis = 1);
    hpa_adjust['Brain'] = average_brain_cols;
    
    return hpa_adjust;

'''
This function imports HPA data (from RNA, IHC, and consensus).
'''

def import_hpa_data():
    hpa_ihc_genes = np.ravel(list(pd.read_csv('../data/HPA/ihc_gene_names.csv', header = None).values))

    hpa_ihc = pd.read_csv('../data/HPA/hpa_ihc_reformat.csv')


    hpa_ihc = hpa_ihc.drop(hpa_ihc.columns[:2], axis= 1)
    hpa_ihc = hpa_ihc.drop('nan', axis= 1)
    hpa_ihc = adjust_hpa_data(hpa_ihc);

    hpa_rna_genes = np.ravel(list(pd.read_csv('../data/HPA/rna_gene_names.csv', header = None).values))

    hpa_rna = pd.read_csv('../data/HPA/hpa_rna_reformat.csv')
    hpa_rna = hpa_rna.drop(hpa_rna.columns[0], axis= 1)
    hpa_rna = adjust_hpa_data(hpa_rna);

    hpa_consensus_full = pd.read_csv('../data/HPA/rna_tissue_consensus.csv');
    hpa_consensus_genes = list(set(hpa_consensus_full['Gene name']))
    brain_parts = ['gangla', 'amygdala', 'pons', 'hippo', 'spinal', 'pituitary', 'cortex', 'thalamus', 'medulla', 'choroid', 'bulb', 'brain', 'white', 'cere']
    hpa_consensus = pd.read_csv('../data/HPA/hpa_consensus_reformat.csv')
    hpa_consensus = hpa_consensus.drop(hpa_consensus.columns[0], axis= 1)
    hpa_consensus = adjust_hpa_data(hpa_consensus, brain_parts)
    
    return hpa_ihc, hpa_ihc_genes, hpa_rna, hpa_rna_genes, hpa_consensus, hpa_consensus_genes


'''
This function imports proteomics data of IPSC-derived EVs from different cell-types.
'''

def import_IPSC_data(brain_celltype, g2u, oldg2u):
    if brain_celltype == None:
        return []
    IPSC_protein_file = '../data/IPSCdata/IPS_'+ brain_celltype + '.csv';
    IPSC_proteins = list(pd.read_csv(IPSC_protein_file).values[:,0])
    IPSC_uniprots = [];

    for protein in IPSC_proteins:
        if protein in g2u.keys():
            IPSC_uniprots.append(g2u[protein]);
            


    return IPSC_uniprots