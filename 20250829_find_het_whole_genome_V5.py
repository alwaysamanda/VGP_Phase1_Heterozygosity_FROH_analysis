#### Script to calculate heterozgyosity based on variants in vcf file ####
## Date: 20250829 (August 29th, 2025)
## Author: Amanda Gardiner
## Version 5 (V3 is 20250728_find_het_calc_V4.py)
## NOTES: This is one part of finding heterozygosity code -- split the original script into two different scripts
##        Did to make code usable within Snakemake snakefile
## NOTES: Updated to include updated path for reading in heterozygosity files from ../groups/.... directory
## NOTES: Updated to read in autosomal chromosome names instead of just counting number of autosomal chromsomes
## NOTES: Updated to exclude all chromosomes <100kb long
## GOAL: Take input of txt file containing chromosomes and variant positions and use it to calculate heterozygosity per chromosome and over whole genome
####

###############################################################################

#### ---- Import necessary libraries ---- ####
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

#### ---- Load in variables from shell ---- ####
num_aut_chr = int(sys.argv[1]) # Number of autosomal chromosomes in the genome
clade = sys.argv[2] # Clade name
full_spec_name = sys.argv[3]
output_per_chr = sys.argv[4]
output_whole_genome = sys.argv[5]
chrom_length_file = sys.argv[6]
chroms = sys.argv[7:(num_aut_chr+7)]

#### ---- FUNCTIONS ---- ####
def read_chromosome_file(chromosome_file):
    with open(chromosome_file, 'r') as file:
        dat = [line.split()[:] for line in file]
    return pd.DataFrame(dat, columns=["chr_name", "chr_length"])

def calc_het(clade, full_species_name, output_file_per_chr, output_file_whole_genome, chroms, chrom_length_file):
    chr_df = read_chromosome_file(chrom_length_file)
    chr_df.iloc[:, 1] = chr_df.iloc[:, 1].astype(int)
    big_chroms = chr_df.loc[(chr_df.iloc[:, 1] > 100000) & (chr_df.iloc[:, 0].isin(chroms)), chr_df.columns[0]].tolist()

    ## Create a dataframe to track mean results for each chromosome
    df = pd.DataFrame({"Chromosomes":big_chroms})
    df['mean_het'] = ''
    df['mean_het_excl_ROH'] = ''

    all_het = []
    all_auto_het = []

    all_het_excl_ROH = []
    all_auto_het_excl_ROH = []

    for chr_name in big_chroms: ## Read in chromosome files
        file_path = os.path.join("../groups", clade, full_species_name, "heterozygosity/Het", chr_name + '_het.txt')
        chrom_file = pd.read_csv(file_path, sep=',', header=0)

        ## Calculate mean het per kb for the chromosome and store it in df
        df.loc[df['Chromosomes']==chr_name,'mean_het'] = np.median(chrom_file['Het_Per_Kb'])
        df.loc[df['Chromosomes']==chr_name,'mean_het_excl_ROH'] = np.median(chrom_file['Het_Per_Kb_excl_ROH'])

        ## Calculate mean het per kb window for the whole genome and just autosomal genome
        all_het.extend(chrom_file['Het_Per_Kb'].tolist())
        all_het_excl_ROH.extend(chrom_file['Het_Per_Kb_excl_ROH'].tolist())
    ## Finish for loop
        ## Save file containing mean heterozygosty per kb window per chromosome
    df.to_csv(output_file_per_chr, header=True, index=False)

    total_mean = np.median(all_het)
    total_mean_excl_ROH = np.median(all_het_excl_ROH)

    for chr_name in big_chroms:
        file_path = os.path.join("../groups", clade, full_species_name, "heterozygosity/Het", chr_name + '_het.txt')
        chrom_file = pd.read_csv(file_path, sep=',')
        all_auto_het.extend(chrom_file['Het_Per_Kb'].tolist())
        all_auto_het_excl_ROH.extend(chrom_file['Het_Per_Kb_excl_ROH'].tolist())
    ## Finish for loop
    total_auto_mean = np.median(all_auto_het)
    total_auto_mean_excl_ROH = np.median(all_auto_het_excl_ROH)

    ## Save results
    results = [
        'The mean heterozygosity per 1kb over the whole genome: \n', total_mean, '\n', 
        'The mean heterozygosity per 1kb over the whole autosomal genome: \n', total_auto_mean, '\n', 
        'The mean heterozygosity per 1kb over the whole genome excluding ROH: \n', total_mean_excl_ROH, '\n',
        'The mean heterozygosity per 1kb over the whole autosomal genome excluding ROH: \n', total_auto_mean_excl_ROH 
    ]

    result_file = open(output_file_whole_genome, 'w')
    result_file.writelines(str(element) for element in results)
    result_file.close()
## End function

run_het_calculations = calc_het(clade, full_spec_name, output_per_chr, output_whole_genome, chroms, chrom_length_file)
