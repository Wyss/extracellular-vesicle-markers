# Finding Candidate Markers for Cell-Type Specific Exosome Isolation From the Brain

## Siddharth Iyer

This repository contains our code to find candidate markers for cell-type specific exosome isolation from the brain. 

### Methods

Broadly, the code uses protein annotations, gene expression data, and proteomics data. We isolate candidate proteins or genes which have transmembrane domains, are specifically expressed in a desired cell-type in the brain, and are specifically expressed in the brain compared to other tissues. Then, we compare to our own mass-spectrometry data of CaptoCore-isolated exosomes from plasma, CSF, and iPS neuron culture, as well as past datasets from the human plasma proteome project (HPPP) and others.

### Data

The data can be found at https://drive.google.com/file/d/1UszDdikP1j9NT8kU1qZmKsViOnuZB3BN/view?usp=sharing. Download this to the same folder where you clone the repository and unzip it there.
Make a new folder called 'results' as well. 

The structure of the repository when set up should be as follows:

neuron_EV_markers/

    ├── bin/           # Scripts and executables

    ├── data/          # Input data files

    ├── figures/       # Figures for the manuscript

    ├── results/       # Analysis outputs

    ├── supp_tables/   # Supplementary tables

./bin/ contains the code for the analysis, ./data/ contains the raw data used in the analysis, ./figures/ contains the figures produced for the paper, ./results contains other figures and results of the code which we found useful to save, and ./supp_tables/ contains the supplementary tables produced for the paper.

### Dependencies

The major Python package requirements and their tested versions are in requirements.txt.

This code was run with Python version 3.7.

### Instructions

To replicate our results, please clone this repository, download the data from the Google Drive link and install the python packages denoted in the requirements.txt file. Then, use the ./bin/figure_producer.ipynb notebook and click "Run All".

If you would like to play around with our code or generate your own lists of candidates we have made a notebook called ./bin/sandbox.ipynb. 

This should output figures in the ./results folder. The figures which were specifically used in the paper should be in the ./figures folder, and supplemental tables should be in the ./supp_tables folder. 

Please let us know if you have difficulty using this!

Siddharth Iyer

iyers@mit.edu

### Questions

For questions, please use the GitHub Discussions forum. For bugs or other problems, please file an issue.
