#### Script to calculate FROH for HGDP analyses ####
## Date: 20260310 (March 10th, 2026)
## Author: Amanda Gardiner
## Version 1
## GOALS: Calculate the FROH normalized by percent of aligned chromosome for HGDP data and append file
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
roh_file = sys.argv[1]
chrom_length_file = sys.argv[2]
output_file_name = sys.argv[3]


#### ---- FUNCTIONS ---- ####
## Function to calculate FROH
def calc_FROH_aut(roh_file, chrom_length_file, output_file_name):
      chrom_dat = pd.read_csv(chrom_length_file, header=None, sep=' ')
      chrom_dat.columns = ['Chrom_name', 'Chrom_length']
      chrom_dat['Chrom_name'] = chrom_dat['Chrom_name'].str.replace('chr', '')
      roh_dat = pd.read_csv(roh_file, sep='\s+')
      roh_dat['length'] = roh_dat['KB']*1000
      chrom_lengths = chrom_dat.set_index('Chrom_name')['Chrom_length']
      roh_dat['FROH'] = roh_dat.iloc[:, 13] / roh_dat.iloc[:, 3].astype(str).map(chrom_lengths)
      roh_dat['FROH_percent'] = roh_dat['FROH']*100

      roh_dat.to_csv(output_file_name, sep=',')
## End function

#### ---- Run function ---- ####
calc_FROH_aut(roh_file, chrom_length_file, output_file_name)
