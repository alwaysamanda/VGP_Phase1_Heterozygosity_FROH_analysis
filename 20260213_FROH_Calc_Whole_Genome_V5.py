#### Script to calculate FROH across the whole genome (Inbreeding Coefficient from ROH Durbin Calculations ####
## Date: 20260213 (February 13th, 2026)
## Author: Amanda Gardiner
## Version 4 (Version 4 is 20250926FROH_Calc_Whole_Genome_V4.py)
## NOTES: Based on rds/users/ag2427/hpc-work/ROH_analyses/20241211_FROH_ROH_Durbin_Calc.R
## NOTES: Doing this to calculate it based on chromosome length from only aligned bases 
## NOTES: Updating this to include ultralong ROH >10Mb
## NOTES: Modified on 20251016 to save the output of aligned size of genome as a separate file
## NOTES: Modified on 20260128 to save the aligned size of each chromosome in the genome size ouput file
## NOTES: Updated to V5 on 20260213 to add in FROH using aln_lengths and also in percent
## NOTES: Updated on 20260320 to recalculate NROH and NROH_aln
####

#### ---- Import necessary libraries ---- ####
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

#### ----- Load in variables from the shell ---- ####
Aln_file = sys.argv[1]
num_aut_chr = int(sys.argv[2])
chromosome_list = sys.argv[3]
roh_data_file = sys.argv[4]
output_file_name = sys.argv[5]
output_align_file = sys.argv[6]
auto_chrom_names = sys.argv[7:(num_aut_chr+7)]

#### ---- FUNCTIONS ---- ####
## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[0] == chromosome]
    return dat

## Define function to calculate sum of the alignments in bases
def get_aln_sum(Alignment_list):
    total_sum = sum(int(line[2]) - int(line[1]) for line in Alignment_list)
    return total_sum
## Finish function

## Function to read chromosome list
def read_chromosome_dat(chromosome_list_file, chromosome_names_list):
    Chrom = pd.read_csv(chromosome_list_file, header=None, sep='\t')
    Chrom.columns = ['Chrom_name', 'Chrom_end']
    if isinstance(chromosome_names_list, str):
            chromosome_names_list = [chromosome_names_list]
    Chrom_autosomal = Chrom[Chrom['Chrom_name'].isin(chromosome_names_list)]
    Chrom_total_length = Chrom['Chrom_end'].sum()
    return Chrom_autosomal, Chrom_total_length
## Finish function

## Function to get sums of ROH in different categories
def get_roh_sum(roh_file_name):
    try:
        dat = pd.read_csv(roh_file_name, header=None)
        dat.columns = ['index', 
                       'chrom', 
                       'start', 
                       'end', 
                       'length', 
                       'aln_length', 
                       'ROH_type', 
                       'ROH_FROH_percent', 
                       'ROH_type_aln', 
                       'ROH_FROH_percent_aln'
                       ]
        dat['length'] = pd.to_numeric(dat['length'], errors='coerce')
        dat['aln_length'] = pd.to_numeric(dat['aln_length'], errors='coerce')
        dat['ROH_FROH_percent'] = pd.to_numeric(dat['ROH_FROH_percent'], errors='coerce')
        dat['ROH_FROH_percent_aln'] = pd.to_numeric(dat['ROH_FROH_percent_aln'], errors='coerce')
        sroh_all_lengths = dat['length'].sum()
        sroh_all_lengths_aln = dat['aln_length'].sum()
        sroh_med_long = dat[dat['length'] >= 500000]['length'].sum()  ## Sum of medium and long ROH lengths
        sroh_med_long_aln = dat[dat['aln_length'] >= 500000]['aln_length'].sum() 
        sroh_long = dat[dat['length'] >= 1000000]['length'].sum() ## Sum of long ROH lengths
        sroh_long_aln = dat[dat['aln_length'] >= 1000000]['aln_length'].sum()
        sroh_ultralong = dat[dat['length'] >= 10000000]['length'].sum() ## Sum of ultralong ROH lengths
        sroh_ultralong_aln = dat[dat['aln_length'] >= 10000000]['aln_length'].sum()

        sroh_one_percent = dat[dat["ROH_FROH_percent"] >= 1]["length"].sum()
        sroh_one_percent_aln = dat[dat["ROH_FROH_percent_aln"] >= 1]["aln_length"].sum()
        sroh_ten_percent = dat[dat["ROH_FROH_percent"] >= 10]["length"].sum()
        sroh_ten_percent_aln = dat[dat["ROH_FROH_percent_aln"] >= 10]["aln_length"].sum()
        sroh_25_percent = dat[dat["ROH_FROH_percent"] >= 25]["length"].sum()
        sroh_25_percent_aln =  dat[dat["ROH_FROH_percent_aln"] >= 25]["aln_length"].sum()
        sroh_50_percent = dat[dat["ROH_FROH_percent"] >= 50]["length"].sum()
        sroh_50_percent_aln = dat[dat["ROH_FROH_percent_aln"] >= 50]["aln_length"].sum()

        nroh_all_lengths = 0 if len(dat) <= 1 else len(dat) - 1 ## Remove one to remove the header from dat
        nroh_all_lengths_aln = 0 if len(dat) <= 1 else ((dat['aln_length'] != 0).sum() - 1) ## Remove the header and only count where aln_length is non-zero
        nroh_med_long = len(dat[dat['length'] >= 500000])
        nroh_med_long_aln = len(dat[dat['aln_length'] >= 500000])
        nroh_long = len(dat[dat['length'] >= 1000000])
        nroh_long_aln = len(dat[dat['aln_length'] >= 1000000])
        nroh_ultralong = len(dat[dat['length'] >= 10000000])
        nroh_ultralong_aln = len(dat[dat['aln_length'] >= 10000000])

        nroh_one_percent = len(dat[dat["ROH_FROH_percent"] >= 1])
        nroh_one_percent_aln = len(dat[dat["ROH_FROH_percent_aln"] >= 1])
        nroh_ten_percent = len(dat[dat["ROH_FROH_percent"] >= 10])
        nroh_ten_percent_aln = len(dat[dat["ROH_FROH_percent_aln"] >= 10])
        nroh_25_percent = len(dat[dat["ROH_FROH_percent"] >= 25])
        nroh_25_percent_aln = len(dat[dat["ROH_FROH_percent_aln"] >= 25])
        nroh_50_percent = len(dat[dat["ROH_FROH_percent"] >= 50])
        nroh_50_percent_aln = len(dat[dat["ROH_FROH_percent_aln"] >= 50])


    except pd.errors.EmptyDataError: ## If there are no ROH found
        print('Note: ROH csv was empty. Skipping.')
        dat = [0]
        sroh_all_lengths = sroh_med_long = sroh_long = sroh_ultralong = 0
        sroh_all_lengths_aln = sroh_med_long_aln = sroh_long_aln = sroh_ultralong_aln = 0
        sroh_one_percent = sroh_ten_percent = sroh_25_percent = sroh_50_percent = 0
        sroh_one_percent_aln = sroh_ten_percent_aln = sroh_25_percent_aln = sroh_50_percent_aln = 0
        nroh_all_lengths = nroh_all_lengths_aln = nroh_med_long = nroh_med_long_aln = nroh_long = nroh_long_aln = nroh_ultralong = nroh_ultralong_aln = 0
        nroh_one_percent = nroh_one_percent_aln = nroh_ten_percent = nroh_ten_percent_aln = nroh_25_percent = nroh_25_percent_aln = nroh_50_percent = nroh_50_percent_aln = 0
    all_roh_sums = pd.DataFrame({
        'ROH_type':[
            'all', 'all_aln', 
            '500kb+', '500kb+_aln', 
            '1Mb+', '1Mb+_aln', 
            '10Mb+', '10Mb+_aln', 
            '1%+', '1%+_aln', 
            '10%+', '10%+_aln', 
            '25%+', '25%+_aln', 
            '50%+', '50%+_aln'
        ], 
        'SROH':[
            sroh_all_lengths, sroh_all_lengths_aln, 
            sroh_med_long, sroh_med_long_aln, 
            sroh_long, sroh_long_aln,
            sroh_ultralong, sroh_ultralong_aln, 
            sroh_one_percent, sroh_one_percent_aln, 
            sroh_ten_percent, sroh_ten_percent_aln, 
            sroh_25_percent, sroh_25_percent_aln, 
            sroh_50_percent, sroh_50_percent_aln
        ], 
        'NROH': [
            nroh_all_lengths, nroh_all_lengths_aln, 
            nroh_med_long, nroh_med_long_aln, 
            nroh_long, nroh_long_aln,
            nroh_ultralong, nroh_ultralong_aln, 
            nroh_one_percent, nroh_one_percent_aln, 
            nroh_ten_percent, nroh_ten_percent_aln, 
            nroh_25_percent, nroh_25_percent_aln, 
            nroh_50_percent, nroh_50_percent_aln
        ]
    })
    
    return all_roh_sums
## Finish function

## Define function to calculate FROH
def calc_FROH_aut(chromosome_list, roh_data_file, Aln_file, output_file_name, output_align_file, auto_chrom_names):
    chrom_dat, total_genome_length = read_chromosome_dat(chromosome_list, auto_chrom_names)
    roh_info = get_roh_sum(roh_data_file)

    all_aln_length = [
        get_aln_sum(read_aln_file(Aln_file, chrom))
        for chrom in chrom_dat.iloc[:, 0]
    ]

    sum_aln_length = sum(all_aln_length)

    with open(output_align_file, "wt") as file:
        proportion_aligned = sum_aln_length / total_genome_length
        file.write(f"The size of the aligned genome is: \n {sum_aln_length} \n")
        file.write(f"The proportion of genome aligned is: \n {proportion_aligned} \n")

        for idx, chrom in enumerate(chrom_dat.iloc[:, 0]):
            file.write(f"The aligned size of {chrom} is: {all_aln_length[idx]}\n")

    cond_aln = roh_info.ROH_type.str.contains("aln") 
    cond_percent = roh_info.ROH_type.str.contains("%")

    conditions = [ 
        cond_aln & ~cond_percent, 
        cond_aln & cond_percent,
        ~cond_aln & cond_percent,  
        ~cond_aln & ~cond_percent 
    ]

    values = [ 
        roh_info['SROH'] / sum_aln_length, 
        roh_info['SROH'] / sum_aln_length, 
        roh_info['SROH'] / total_genome_length, 
        roh_info['SROH'] / total_genome_length 
    ]

    roh_info['FROH'] = np.select(conditions, values)
    roh_info['FROH_percent'] = roh_info['FROH'] * 100

    roh_info.to_csv(output_file_name, sep=",")
## Finish Function

#### ---- Run Function ---- ####
calc_FROH_aut(chromosome_list, roh_data_file, Aln_file, output_file_name, output_align_file, auto_chrom_names)

