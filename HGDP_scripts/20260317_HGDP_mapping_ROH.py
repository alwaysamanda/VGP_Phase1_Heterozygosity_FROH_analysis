#### Script to convert ROH from bases to cM ####
## Date: 20260317 (March 17th, 2026)
## Version: 1
## Author: Amanda Gardiner
## GOAL: Use genetic map to convert HGDP ROH positions from bp to cM in size
####

#### ---- Import necessary libraries ---- ####
import sys
from datetime import date
import pandas as pd
import numpy as np

#### ---- Load in data and variables ---- ####
args = sys.argv
date = date.today()
data = pd.read_csv(args[1])
genetic_map = pd.read_csv(
    args[2], 
    sep=' ')
output_file_name = args[3]

#### ---- Sort the genetic map ---- ####
genetic_map = genetic_map.sort_values(['chr', 'position']).reset_index(drop=True)

#### ---- Interpolate bp to cM for start positions of ROH ---- ####
data['POS1_cm'] = np.nan
data['POS2_cm'] = np.nan

for chrom, chrom_map in genetic_map.groupby('chr'):
    mask = data['CHR'] == chrom
    if mask.sum() == 0:
        continue
    data.loc[mask, 'POS1_cm'] = np.interp(
        data.loc[mask, 'POS1'],
        chrom_map['position'],
        chrom_map['Genetic_Map(cM)']
    )
    data.loc[mask, 'POS2_cm'] = np.interp(
        data.loc[mask, 'POS2'],
        chrom_map['position'],
        chrom_map['Genetic_Map(cM)']
    )

data['length_cm'] = data['POS2_cm'] - data['POS1_cm']

#### ---- Save output ---- ####
data.to_csv(output_file_name, index=False, header=True)



