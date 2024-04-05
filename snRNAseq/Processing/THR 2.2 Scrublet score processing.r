# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created.

# 2.2 Scrublet score processing - calculate Scrublet scores for each nuclei and save the output in folders per sample

# Author: Thomas Rust (THR) adapted from Mirjam Koster (MKO) and Astrid Alsema 
# Date 26/10/2023 (updated)

# Run on Jupyter Lab in a Python 3 ipykernel

# Load python packages
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas
import glob
from platform import python_version

# set output directory
print("Current Python Version-", python_version())
path = os.getcwd() +"/Output/"
print("Current path is ", path)

# Run Scrublet for each sample in the object list and save the output
samplename = []
for dir in glob.glob(path + "R_QC/Scrublet/counts/*/matrix.mtx"):
    fileName = dir.split('/')[10] # change the number
    print(' print the sample id:')
    print(fileName) # should print the sample , otherwise change the number
    samplename.append(fileName)

for i in range(0, len(samplename)):
    input_dir =  path + "R_QC/Scrublet/counts/" + samplename[i]
    counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
    genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t'))
    print('running scrublet...')
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.10)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
    np.savetxt(path + "R_QC/Scrublet/" + samplename[i]  + '_doublet_scores.csv', doublet_scores)
    np.savetxt(path + "R_QC/Scrublet/" + samplename[i]  + '_predicted_doublets.csv', predicted_doublets)
    
# We now have the doublet scores for both DoubletFinder and Scrublet. Compare the output (particularly the UMAPs for identification of doublets 
# mainly between astrocytes and microglia. We decided to use Scrublet to remove Doublets from the dataset. 
# Proceed to THR 2.3 Doublet exclusion removal