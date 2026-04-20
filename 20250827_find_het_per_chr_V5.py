#### Script to calculate heterozgyosity based on variants in vcf file ####
## Date: 20250827 (June 27th, 2025)
## Author: Amanda Gardiner
## Version 5 (V4 is 20250618_find_het_calc_V4.R)
## GOAL: Take input of txt file containing chromosomes and variant positions and use it to calculate heterozygosity per chromosome and over whole genome
## NOTES: Modifying this to make the alternate window size 100kb and to exclude chromosomes smaller than 100kb from analysis
## NOTES: Modifying this to use the individual chrom var files instead of the larger ones
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
from cyvcf2 import VCF

## Load in variables from shell
chromosome_length_file = sys.argv[1] # File containing chromosome names and their lengths
clade = sys.argv[2] # Clade name
full_species_name = sys.argv[3] #Species name/reference genome name -- directory for storing results
window_length = int(sys.argv[4]) # Window length defined for calculating het
window_interval = int(sys.argv[5]) # How often windows start, defines their level of overlap
roh_file = sys.argv[6] # File containing ROH start, end, and length for whole genome
chrom_name = sys.argv[7]
output = sys.argv[8]
var_file = VCF(sys.argv[9])


#### ---- FUNCTIONS ---- ####
def filter_chr_file_by_chrom(chromosome_file, chromosome):
    with open(chromosome_file, 'r') as file:
        dat = [word for line in file if line.split()[0] == chromosome for word in line.split()]
    return dat

def filter_roh_by_chrom(roh_data, chromosome):
    dat=[]
    with open(roh_data, 'r') as file:
        for line in file:
            splitline = line.split(',')
            if splitline[0] == chromosome:
                splitline[3] = splitline[3].replace('\n', '')
                dat.append(splitline)
    return pd.DataFrame(dat)

def read_vcf_file(vcf_file, chrom_name):
    target_var_pos = []
    for variant in vcf_file(chrom_name):
        target_var_pos.append(variant.POS)
    return np.array(target_var_pos)

def calc_het(roh_file, chromosome_length_file, window_length, window_interval, clade, full_species_name, chrom, output, var_file):
    chrom_dat = filter_chr_file_by_chrom(chromosome_length_file, chrom)
    chrom_name = chrom_dat[0]
    chrom_length = int(chrom_dat[1])
    if chrom_length < 100000:
        ## This chromosome will be excluded from analysis due to small size
        with open(output, "wt") as file:
            file.write("NA")
            file.write(f"{chrom_name} is less than 100kb and too short to analyze")
            file.close()
    else:
        if chrom_length < window_length:
            window_length_alt = 100000
            window_interval_alt = window_length_alt/2
            window_starts = np.arange(0, (chrom_length-window_length_alt+1), window_interval_alt)
            window_ends = window_starts + window_length_alt
            window_sizes = np.repeat(window_length_alt, len(window_starts))
        else:
            window_starts = np.arange(0, (chrom_length-window_length+1), window_interval)
            window_ends = window_starts + window_length 
            window_sizes = np.repeat(window_length, len(window_starts))
        # single_chr_df = pd.read_csv(var_file, sep="\t", header=None)
        variant_positions = read_vcf_file(var_file, chrom_name)
        het = np.zeros(len(window_starts)) ## Create heterozygosity vector
        het_wo_roh = np.zeros(len(window_starts)) ## heterozygosity vector excluding ROH
        single_chr_results = pd.DataFrame({
            'Start': window_starts, 
            'End': window_ends, 
            'Het': het, 
            'Het_excl_ROH': het_wo_roh, 
            'Window_Size': window_sizes, 
            'Window_Size_excl_ROH': window_sizes
        }) ## Create dataframe to store results for the given chromosome
        ## Count the number of variants in each window
        starts = single_chr_results.iloc[:, 0].astype(int)
        ends = single_chr_results.iloc[:, 1].astype(int)
        # variant_positions = single_chr_df.iloc[:, 2].astype(int)
        single_chr_results.loc[:, 'Het'] = [
            ((variant_positions >= start) & (variant_positions <= end - 1)).sum()
            for start, end in zip(starts, ends)
            ]
        # single_chr_results.loc[:,'Het'] = [(variant_positions.between(start, end - 1)).sum() for start, end in zip(starts, ends)] ## Sum all variants in each window, not worrying about ROH
        ## Count number of variants in each window, excluding those found in ROH
        ## Get ROH positions
        roh_dat = filter_roh_by_chrom(roh_file, chrom_name)
        ## If there are ROH present on chromosome -- Calculate het
        if roh_dat.empty == False:
            ROH_starts = (roh_dat.iloc[:, 1].values).astype(int) 
            ROH_ends = (roh_dat.iloc[:, 2].values).astype(int)    
            # variant_positions = (single_chr_df.iloc[:, 2].values).astype(int)  
            for k in range(len(single_chr_results)):  
                start, end = single_chr_results.iloc[k, [0, 1]].astype(int) ## Start and end of a given window
                var_in_win = variant_positions[(variant_positions >= start) & (variant_positions < end)] ## Find all variants within a given window
                ## Find the first ROH that starts after 'start'
                # l = int(np.searchsorted(ROH_starts, start, side='right') - 1)
                l = np.where(ROH_starts > start)[0]
                if l.size > 0:
                    l = int(l[0])
                else:
                    l = -1
                if l < 0 or end < ROH_starts[l]:  
                    ## No ROH overlap, count all variants in the range
                    sum_variants = len(var_in_win)
                    single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
                elif start >= ROH_starts[l] and end <= ROH_ends[l]:  
                    ## Window fully inside ROH
                    sum_variants = 0  
                    single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
                elif start < ROH_starts[l] and end > ROH_starts[l] and end <= ROH_ends[l]:  
                    ## Window overlaps start of ROH
                    var_in_win_excl_ROH = var_in_win[(var_in_win < ROH_starts[l])]
                    sum_variants = len(var_in_win_excl_ROH)
                    single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants  
                    single_chr_results.loc[k,'Window_Size_excl_ROH'] = ROH_starts[l] - start
                elif start >= ROH_starts[l] and start <= ROH_ends[l] and end > ROH_ends[l]:  
                    ## Window overlaps end of ROH
                    var_in_win_excl_ROH = var_in_win[(var_in_win > ROH_ends[l])]
                    sum_variants = len(var_in_win_excl_ROH)
                    single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
                    single_chr_results.loc[k,'Window_Size_excl_ROH'] = end - ROH_ends[l]
                elif start < ROH_starts[l] and end > ROH_ends[l]:  
                    ## Window Completely contains ROH
                    variants_pre_ROH = var_in_win[(var_in_win >= start) & (var_in_win < ROH_starts[l])]
                    sum_variants_pre_ROH = len(variants_pre_ROH)
                    variants_post_ROH = var_in_win[(var_in_win > ROH_ends[l]) & (var_in_win < end)]
                    sum_variants_post_ROH = len(variants_post_ROH)
                    sum_variants = sum_variants_pre_ROH + sum_variants_post_ROH
                    single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
                    single_chr_results.loc[k,'Window_Size_excl_ROH'] = (ROH_starts[l] - start) + (end - ROH_ends[l])
            ## Finish for loop
        elif roh_dat.empty == True: ## Calculate het if there are NO ROH on chr
            single_chr_results.loc[:,'Het_excl_ROH'] = single_chr_results.loc[:,'Het'] ## Het_excl_ROH is same as het given that there are no ROH on chromosome
        ## Calculate heterozygosity per kb
        single_chr_results['Het_Per_Kb'] = (single_chr_results['Het']/single_chr_results['Window_Size'])*1000
        single_chr_results['Het_Per_Kb_excl_ROH'] = (single_chr_results['Het_excl_ROH']/single_chr_results['Window_Size_excl_ROH'])*1000
        ## Add chromosome column
        chrom_list=[chrom_name]*single_chr_results.shape[0]
        single_chr_results.insert(loc = 0, column = 'chrom', value=chrom_list)
        ## Save results for chromosome
        single_chr_results.to_csv(output, index=False, header=True)
    ## End function

run_het_calculations = calc_het(roh_file, chromosome_length_file, window_length, window_interval, clade, full_species_name, chrom_name, output, var_file)
