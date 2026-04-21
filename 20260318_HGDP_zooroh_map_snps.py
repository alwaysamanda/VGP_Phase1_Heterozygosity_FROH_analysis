#### Script to interpolate HGDP SNPs in centimorgans ####
## Date: 20260318 (March 18th, 2026)
## Version: 1
## Author: Amanda Gardiner
## GOAL: Get position of HGDP SNPs in centimorgans for running in RZooROH
####

#### ---- Import necessary libraries ---- ####
import sys
from datetime import date
import pandas as pd
import numpy as np

#### ---- Load in data and variables ---- ####
args = sys.argv
chrom = int(args[1])
gmap = pd.read_csv(args[2], sep=' ')
snps = pd.read_csv(args[3], sep='\t', names=['snp_id','chr','pos_bp'])
output_file = args[4]

#### ---- Interpolate position ---- ####
gmap_chr = gmap[gmap['chr'] == chrom]

snps['pos_cm'] = np.interp(
    snps['pos_bp'],
    gmap_chr['position'],
    gmap_chr['Genetic_Map(cM)'])

#### ---- Save file ---- ####
snps[['snp_id','chr','pos_cm']].to_csv(output_file, sep='\t', index=False, header=False)