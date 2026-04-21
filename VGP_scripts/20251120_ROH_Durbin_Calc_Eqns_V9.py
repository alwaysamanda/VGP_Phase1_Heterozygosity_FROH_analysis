#### WRITING DURBIN'S EQUATIONS TO CALCULATE ROH FOR VGP GENOMES ####
## Date: 20251120 (November 20th, 2025)
## Author: Amanda Gardiner
## Version 9 (V8 is 20251001_ROH_Durbin_Calc_Eqns_V8.py)
## NOTES: 
##          Start with X0 = 0, a0 = b0 = c0 = S0 = 0
##          ai = start, bi = end, ci = difference
##          Move along chromosome between SNPs
##          Si = Si-1 + (Xi - Xi-1) - P
##          Si-1 = score before going to this region
##          Xi - Xi-1 = distance between this SNP and previous SNP
##          P = set penalty (try 1,000,000 and then 100,000)
##          If Si ≤ 0, then set ai = Xi, bi = Xi, ci = 0
##          If Si > 0, then if Si > ci, make ci = Si, and bi = Xi
##          Repeat along whole chromosome to get regions of interest

## NOTES: Doing this to modify input from ALN file to account for swapped position
##        of Ref Chr names in Aln files using Richard Durbin's FastGA alingments of VGP data freeze
## NOTES: Updating this to avoid error with last version where ROH could overlap -- which is not possible
## NOTES: Updating this to take a vcf file as input for the SNVs instead of a txt file
####

#### ---- Import necessary libraries ---- ####
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os
from cyvcf2 import VCF

#### ----  Load in variables from the shell script ---- ####
chrom = sys.argv[1] ## Chromosome ROH is being calculated for
Aln_file = sys.argv[2] ## File containing alignment of alternate to reference
Var_file = VCF(sys.argv[3]) ## File containing variants
output_name = sys.argv[4] ## Name of output file

#### ----  FUNCTIONS ---- ####
## Function to read variants for a given chromosome and store it in an array
def read_vcf_file(vcf_file, chromosome):
    target_var_pos = []
    for variant in vcf_file(chromosome):
        target_var_pos.append(variant.POS)
    return np.array(target_var_pos)

## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[0] == chromosome]
    return dat

## Function to read alignment bed file
def aln_map_intervals(Alignment_list):
    starts = np.array([int(line[1]) for line in Alignment_list], dtype=np.int64) ## Start pos of each aln
    ends   = np.array([int(line[2]) for line in Alignment_list], dtype=np.int64) ## End pos of each aln
    lengths = ends - starts ## Length of alns
    
    ## Sort intervals by start
    order = np.argsort(starts)
    starts, ends, lengths = starts[order], ends[order], lengths[order]
    
    ## Get cumulative aligned length 
    cum_lengths = np.cumsum(lengths)
    
    return starts, ends, cum_lengths

## Function to merge ROH if they are overlapping
def merge_roh(df_roh):
    df_roh = df_roh.sort_values(by="start").reset_index(drop=True)
    merged = []
    
    for _, row in df_roh.iterrows():
        if not merged:
            merged.append(row)
        else:
            prev = merged[-1]
            # If current ROH overlaps or touches the previous one
            if row['start'] <= prev['end']:
                # Merge them by extending the end
                prev['end'] = max(prev['end'], row['end'])
            else:
                merged.append(row)
    
    return pd.DataFrame(merged)

## Function to get the sum of aligned bases for each variant
def compute_var_aln_sum(var_pos, starts, ends, cum_lengths):
    indices = np.searchsorted(starts, var_pos, side="right") - 1 ## Switch to 0 indexing
    
    var_aln_sum = np.zeros(len(var_pos), dtype=np.int64) ## var_aln_sum == sum of all alingned bases at the variant's point
    for i, idx in enumerate(indices):
        if idx >= 0:
            ## Sum of aligned bases up to end of previous alignment
            prev_sum = cum_lengths[idx-1] if idx > 0 else 0
            
            ## Add partial length of alignment to that point
            if starts[idx] <= var_pos[i] < ends[idx]:
                var_aln_sum[i] = prev_sum + (var_pos[i] - starts[idx] + 1)
            else:
                var_aln_sum[i] = cum_lengths[idx]
        else:
            var_aln_sum[i] = 0
    return var_aln_sum

## Function to calculate positions of Runs of Homozygosity in a given chromosome
def calculate_ROH(chrom, Aln_file, Var_file):
    aln_list = read_aln_file(Aln_file, chrom)
    starts, ends, cum_lengths = aln_map_intervals(aln_list)
    var_pos = read_vcf_file(Var_file, chrom)
    ## Compute alignment sums at variant positions
    var_aln_sum = compute_var_aln_sum(var_pos, starts, ends, cum_lengths)
    ## Build DataFrame
    df = pd.DataFrame({
        'variant': var_pos,
        'var_aln_sum': var_aln_sum,
        'score': np.zeros(len(var_pos), dtype=int)
    })
    # Apply scoring equation
    Penalty = 1e5
    df['score'] = np.maximum(
        0,
        (np.roll(df['score'].astype(int), 1)
        + (df['var_aln_sum'].astype(int) - np.roll(df['var_aln_sum'], 1))
        - Penalty)
    )
    ## Extract the start and end points for each ROH
    roh_boundaries = []
    in_roh = False
    start_idx = None
    for i in range(len(df)):
        if df.loc[i, 'score'] > 0:
            if not in_roh:
                in_roh = True ## Start a new ROH if not already in one
                start_idx = i
        else:
            if in_roh:
                end_idx = i - 1 ## End ROH when score is negative/0
                start_variant = df.loc[start_idx - 1, 'variant'] if start_idx > 0 else df.loc[start_idx, 'variant']
                end_variant = df.loc[end_idx, 'variant']
                roh_boundaries.append((start_variant, end_variant))
                in_roh = False
    if in_roh:
        end_idx = len(df) - 1 # If final ROH goes through the end of the chromosome
        start_variant = df.loc[start_idx - 1, 'variant'] if start_idx > 0 else df.loc[start_idx, 'variant']
        end_variant = df.loc[end_idx, 'variant']
        roh_boundaries.append((start_variant, end_variant))
    # If no ROH found
    if len(roh_boundaries) == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'length'])
    result = pd.Series(roh_boundaries)
    ## Use start and end points to create a data frame containing ROH start, end, and length
    ROH_start = [roh[0] for roh in result]
    ROH_end = [roh[1] for roh in result]
    df_roh = pd.DataFrame({
        'start': ROH_start, 
        'end': ROH_end
    })
    df_roh = merge_roh(df_roh)
    df_roh['length'] = df_roh['end'] - df_roh['start']
    chrom_list = [chrom]*df_roh.shape[0]
    df_roh.insert(loc = 0, column = 'chrom', value=chrom_list)
    return df_roh
#### ---- END ---- ####

#### ---- Find ROH and save them in output csv ---- ####
ROH_report = calculate_ROH(chrom, Aln_file, Var_file)

if ROH_report.empty == True:
    print('No ROH detected in', chrom)

ROH_report.to_csv(output_name, index=False, header=True)







# #### ----  FUNCTIONS ---- ####
# ## Define function to read the chromosome length file and get length
# def read_length_file(chromosome_length_file, chromosome):
#     with open(chromosome_length_file, 'r') as file:
#         dat = [line.split()[1] for line in file if line.split()[0] == chromosome]
#         dat = int(dat[0])
#     return dat

# ## Define function to read the alignment file and store it in an array
# def read_aln_file(Alignment_file, chromosome):
#     with open(Alignment_file, 'r') as file:
#         dat = [line.split() for line in file if line.split()[1] == chromosome]
#     return dat

# ## NOTE -- look into making this function faster/more efficient as paftools should already only have uniquely aligned regions
# ## Define function to parse through whole chromosome and map each individual base and whether it is aligned or not
# def aln_map(Alignment_list, chrom_length):
#     base_sums = np.zeros(chrom_length, dtype=int) ## Create an array the length of the chromosome
#     for start, end in ((int(line[2]), int(line[3])) for line in Alignment_list):
#         base_sums[start:end] += 1 ## Calculate depth of alingments at each individual base
#     base_gross_cov = np.where(base_sums > 0, 1, 0) ## Convert depth to binary of whether it is covered or not
#     aln_sums = np.cumsum(base_gross_cov) ## Get alingment sum at each given base
#     return np.vstack((np.arange(1, chrom_length + 1), base_sums, aln_sums.astype(int)))
# ## Finish function

# def aln_map(Alignment_list, chrom_length):
#     base_sums = np.zeros(chrom_length, dtype=np.uint8)
#     for start, end in ((int(line[2]), int(line[3])) for line in Alignment_list):
#         base_sums[start:end] = 1 ## Calculate depth of alingments at each individual base
#     aln_sums = np.cumsum(base_gross_cov)

# ## Function to read variants for a given chromosome and store it in an array
# def read_vcf_file(vcf_file, chromosome):
#     target_var_pos = []
#     for variant in vcf_file(chromosome):
#         target_var_pos.append(variant.POS)
#     return np.array(target_var_pos)

# ## Function to merge ROH if they are overlapping
# def merge_roh(df_roh):
#     df_roh = df_roh.sort_values(by="start").reset_index(drop=True)
#     merged = []
    
#     for _, row in df_roh.iterrows():
#         if not merged:
#             merged.append(row)
#         else:
#             prev = merged[-1]
#             # If current ROH overlaps or touches the previous one
#             if row['start'] <= prev['end']:
#                 # Merge them by extending the end
#                 prev['end'] = max(prev['end'], row['end'])
#             else:
#                 merged.append(row)
    
#     return pd.DataFrame(merged)

# ## Define function to calculate ROH
# def calculate_ROH(chrom_length_file, chrom, Aln_file, Var_file):
#     chrom_length = read_length_file(chrom_length_file, chrom)
#     aln_list = read_aln_file(Aln_file, chrom)
#     base_aln_map = aln_map(aln_list, chrom_length)
#     var_pos = read_var_file(Var_file, chrom)

#     score = np.zeros(len(var_pos), dtype=int) ## Create array to store the score at each base position
#     aln_sums = np.zeros(len(var_pos), dtype=int) ## Create array to store the total number of aligned bases at each position
#     df = pd.DataFrame({'variant': var_pos, 
#                        'var_aln_sum': aln_sums, 
#                        'score': score})
    
#     indices = np.searchsorted(base_aln_map[0, :], df['variant']) - 1 ## Find the variants in the base_aln_map and subtract one to get them to 0-indexing
#     df['var_aln_sum'] = base_aln_map[2, indices] ## var_aln_sum == sum of all alingned bases at the variant's point

#     ## Run equation to calculate score based on aligned chr length (var_aln_sum)
#     Penalty = 1e+05 ## Define penalty -- currently 100,000
#     df['score'] = np.maximum(0, (np.roll(df['score'].astype(int), 1) + (df['var_aln_sum'].astype(int) - np.roll(df['var_aln_sum'], 1)) - Penalty))

#     ## Extract the start and end poitns for each ROH
#     roh_boundaries = []
#     in_roh = False
#     start_idx = None

#     for i in range(len(df)):
#         if df.loc[i, 'score'] > 0:
#             if not in_roh:
#                 in_roh = True ## Start a new ROH if not already in one
#                 start_idx = i
#         else:
#             if in_roh:
#                 end_idx = i - 1 ## End ROH when score is negative/0
#                 start_variant = df.loc[start_idx - 1, 'variant'] if start_idx > 0 else df.loc[start_idx, 'variant']
#                 end_variant = df.loc[end_idx, 'variant']
#                 roh_boundaries.append((start_variant, end_variant))
#                 in_roh = False

#     if in_roh:
#         end_idx = len(df) - 1 # If final ROH goes through the end of the chromosome
#         start_variant = df.loc[start_idx - 1, 'variant'] if start_idx > 0 else df.loc[start_idx, 'variant']
#         end_variant = df.loc[end_idx, 'variant']
#         roh_boundaries.append((start_variant, end_variant))

#     # If no ROH found
#     if len(roh_boundaries) == 0:
#         return pd.DataFrame(columns=['chrom', 'start', 'end', 'length'])

#     result = pd.Series(roh_boundaries)

#     ## Use start and end points to create a data frame containing ROH start, end, and length
#     ROH_start = [roh[0] for roh in result]
#     ROH_end = [roh[1] for roh in result]
#     df_roh = pd.DataFrame({
#         'start': ROH_start, 
#         'end': ROH_end
#     })

#     df_roh = merge_roh(df_roh)

#     df_roh['length'] = df_roh['end'] - df_roh['start']
#     chrom_list = [chrom]*df_roh.shape[0]
#     df_roh.insert(loc = 0, column = 'chrom', value=chrom_list)

#     return df_roh
# #### ---- END ---- ####

# #### ---- Find runs of homozygosity ---- ####
# ROH_report = calculate_ROH(chrom_length_file, chrom, Aln_file, Var_file)

# if ROH_report.empty == True:
#     print('No ROH detected in', chrom)

# ROH_report.to_csv(output_name, index=False, header=True)

