#### Script to calculate FROH (Inbreeding Coefficient from ROH Durbin Calculations ####
## Date: 20260213 (February 13th, 2026)
## Author: Amanda Gardiner
## Version 5 (Version 4 is 20250925_FROH_per_chr_calc_V4.py)
## GOAL: Calculate The Inbreeding coefficient per chromosome for an individual,
##       and report results both in orignial format and normalized for the size of the chromosome so that they can be compared
## NOTES: Initial script pulled from 20250106_FROH_Calc.R
## NOTES: Each time I am running this script, I am running it for a specific chromosome -- as such I do not need to loop through all of them
## NOTES: Altering this from the previous script to get the array of aligned bases and calculate l_aut and chr length based on aligned bases only
## NOTES: Altering this to calculate the chromosome length within the script rather than passing it as a separate variable
## NOTES: Adding in category for ultralong ROH
## NOTES: Altering from previous version to calculate length with total chrom and aligned length with aligned chrom
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
Roh_file = sys.argv[2]
chrom = sys.argv[3]
chrom_length_file = sys.argv[4]
output_file_name = sys.argv[5]

#### ---- FUNCTIONS ---- ####
## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[0] == chromosome]
    return dat

## Define function to get total length of chromosome
def get_chrom_length(chrom_length_file, chromosome):
    df = pd.read_csv(chrom_length_file, header=None, sep='\t')
    df.columns = ['Chrom_name', 'Chrom_end']
    result = df.loc[df['Chrom_name'] == chromosome, 'Chrom_end']
    
    if result.empty:
        raise ValueError(f"Chromosome '{chromosome}' not found in file")
    
    chrom_length = result.iloc[0]

    return chrom_length

## Define function to calculate sum of the alignments in bases
def get_aln_sum(Alignment_list):
    total_sum = sum(int(line[2]) - int(line[1]) for line in Alignment_list)
    return total_sum
## Finish function

## Define function to read the ROH data
def read_ROH_dat(roh_file_name):
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
    except pd.errors.EmptyDataError:
        print('Note: ROH csv was empty. Skipping.')
        dat = pd.DataFrame(
            columns=[
                'index', 
                'chrom', 
                'start', 
                'end', 
                'length', 
                'aln_length', 
                'ROH_type', 
                'ROH_FROH_percent', 
                'ROH_type_aln', 
                'ROH_FROH_percent_aln'
                ])
    return dat
## Finish function

## Define function to calculate FROH
def calc_FROH(chrom, Aln_file, Roh_file, chrom_length_file, output_file_name):
    aln_list = read_aln_file(Aln_file, chrom)
    total_aln_sum = get_aln_sum(aln_list)
    chrom_length = get_chrom_length(chrom_length_file, chrom)
    roh_data = read_ROH_dat(Roh_file)

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
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0
        ], 
        'NROH':[
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0, 
            0, 0
        ]
    })

    if not roh_data.empty:
        chr_roh = roh_data[roh_data["chrom"] == chrom]
        all_roh_sums['SROH'] = [
            chr_roh['length'].sum(), 
            chr_roh['aln_length'].sum(), 
            chr_roh[chr_roh['length'] >= 500000]['length'].sum(), 
            chr_roh[chr_roh['length'] >= 500000]['aln_length'].sum(), 
            chr_roh[chr_roh['length'] >= 1000000]['length'].sum(), 
            chr_roh[chr_roh['length'] >= 1000000]['aln_length'].sum(), 
            chr_roh[chr_roh['length'] >= 10000000]['length'].sum(), 
            chr_roh[chr_roh['length'] >= 10000000]['aln_length'].sum(), 
            chr_roh[chr_roh["ROH_FROH_percent"] >= 1]["length"].sum(),
            chr_roh[chr_roh["ROH_FROH_percent"] >= 1]["aln_length"].sum(),
            chr_roh[chr_roh["ROH_FROH_percent"] >= 10]["length"].sum(),
            chr_roh[chr_roh["ROH_FROH_percent"] >= 10]["aln_length"].sum(),
            chr_roh[chr_roh["ROH_FROH_percent"] >= 25]["length"].sum(), 
            chr_roh[chr_roh["ROH_FROH_percent"] >= 25]["aln_length"].sum(), 
            chr_roh[chr_roh["ROH_FROH_percent"] >= 50]["length"].sum(),
            chr_roh[chr_roh["ROH_FROH_percent"] >= 50]["aln_length"].sum()
        ]

        all_roh_sums['NROH'] = [
            len(chr_roh), 
            len(chr_roh),
            len(chr_roh[chr_roh['length'] >= 500000]), 
            len(chr_roh[chr_roh['aln_length'] >= 500000]), 
            len(chr_roh[chr_roh['length'] >= 1000000]), 
            len(chr_roh[chr_roh['aln_length'] >= 1000000]), 
            len(chr_roh[chr_roh['length'] >= 10000000]), 
            len(chr_roh[chr_roh['aln_length'] >= 10000000]), 
            len(chr_roh[chr_roh["ROH_FROH_percent"] >= 1]), 
            len(chr_roh[chr_roh["ROH_FROH_percent_aln"] >= 1]), 
            len(chr_roh[chr_roh["ROH_FROH_percent"] >= 10]), 
            len(chr_roh[chr_roh["ROH_FROH_percent_aln"] >= 10]), 
            len(chr_roh[chr_roh["ROH_FROH_percent"] >= 25]), 
            len(chr_roh[chr_roh["ROH_FROH_percent_aln"] >= 25]), 
            len(chr_roh[chr_roh["ROH_FROH_percent"] >= 50]), 
            len(chr_roh[chr_roh["ROH_FROH_percent_aln"] >= 50])
        ]

        cond_aln = all_roh_sums.ROH_type.str.contains("aln") 
        cond_percent = all_roh_sums.ROH_type.str.contains("%")

        conditions = [ 
            cond_aln & ~cond_percent, 
            cond_aln & cond_percent,
            ~cond_aln & cond_percent,  
            ~cond_aln & ~cond_percent 
         ]

        values = [ 
            all_roh_sums['SROH'] / total_aln_sum, 
            all_roh_sums['SROH'] / total_aln_sum, 
            all_roh_sums['SROH'] / chrom_length, 
            all_roh_sums['SROH'] / chrom_length 
        ]

        all_roh_sums['FROH'] = np.select(conditions, values)
        all_roh_sums['FROH_percent'] = all_roh_sums['FROH'] * 100



    all_roh_sums.to_csv(output_file_name,  sep=",")
## End function

#### ---- Run Function ---- ####
calc_FROH(chrom, Aln_file, Roh_file, chrom_length_file, output_file_name)