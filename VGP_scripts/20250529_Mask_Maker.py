#### Script to create mask files for MSMC input to exclude any non-aligned regions ####
## Date: 20250528 (May 29th, 2025)
## Author: Amanda Gardiner
## Version 1
## GOAL: Take alignment of haplotypes and only use regions with alignments, adding them to a bed file
## NOTES: Output file will be in bed format
##        ChromStart in column 2 -- (the first base on the chromosome is numbered 0 i.e. the number is zero-based)
##        ChromEnd in column 3 -- This position is non-inclusive, unlike chromStart (the first base on the chromosome is numbered 1 i.e. the number is one-based).

####

###############################################################################

#### Import necessary libraries ####
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

#### Load in necessary data ####
chrom = sys.argv[1]
chrom_length_file = sys.argv[2]
clade = sys.argv[3]
Aln_file = sys.argv[4]
output_file_name = sys.argv[5]

#### Define Functions ####
## Function to read the txt file with chromosome lengths
def filter_chr_file_by_chrom(chromosome_file, chromosome):
    with open(chromosome_file, 'r') as file:
        dat = [word for line in file if line.split()[0] == chromosome for word in line.split()]
    return dat

## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[1] == chromosome]
    return dat

## Define function to parse through whole chromosome and map each individual base and whether it is aligned or not
def aln_map(Alignment_list, chrom_length):
    base_sums = np.zeros(chrom_length, dtype=int) ## Create an array the length of the chromosome
    for start, end in ((int(line[2]), int(line[3])) for line in Alignment_list):
        base_sums[start:end] += 1 ## Calculate depth of alingments at each individual base
    base_gross_cov = np.where(base_sums > 0, 1, 0) ## Convert depth to binary of whether it is covered or not
    aln_sums = np.cumsum(base_gross_cov) ## Get alingment sum at each given base
    return np.vstack((np.arange(1, chrom_length + 1), base_sums, aln_sums.astype(int)))
## Finish function

## Define function to make the mask for the chromosome
def make_mask(aln_array, chromosome, output_file):
    positions = aln_array[0]
    aln = aln_array[1]
    result = []
    in_range = False
    start = None
    for i in range(len(aln)):
        if aln[i] == 1:
            if not in_range:
                start = pd.to_numeric(positions[i])-1 ## Get 0-indexing to make start inclusive but end not since they are on different indexing
                in_range = True
        else:
            if in_range:
                end = positions[i]
                result.append((start, end))
                in_range = False
    if in_range:
        end = positions[-1]
        result.append((start, end))
    ## Write output file
    with open(output_file, 'w') as f:
        for start, end in result:
            f.write(f"{chromosome}\t{start}\t{end}\n")
    print(f"Saved {len(result)} ranges to {output_file}")
## Finish function
#### All Functions Defined####

chr_info = filter_chr_file_by_chrom(chrom_length_file, chrom)
chrom_length = int(chr_info[1])

aln_list = read_aln_file(Aln_file, chrom)
# print(aln_list)

base_aln_map = aln_map(aln_list, chrom_length)

make_mask(base_aln_map, chrom, output_file_name)

