#### Script to get ROH bin population averages for ROH in BP and cM ####
## Date: 20260317 (March 17th, 2026)
## Version: 1
## Author: Amanda Gardiner
## GOAL: Get the average for each population for porportion of roh in each bin with both measures
####

#### ---- Import necessary libraries ---- ####
import sys
from datetime import date
import pandas as pd
import numpy as np

#### ---- Load in data and variables ---- ####
args = sys.argv
bp_bins = pd.read_csv(args[1])
cm_bins = pd.read_csv(args[2])
bp_output = args[3]
cm_output = args[3]

#### ---- Get averages for bins in bp ---- ####
columns=[col for col in bp_bins.columns if col.startswith('bp_log_bin_')]
means_bp_df = bp_bins.groupby('Population_elastic_ID')[columns].mean().reset_index()

#### ---- Get averages for bins in cM ---- ####
columns=[col for col in cm_bins.columns if col.startswith('cm_log_bin_')]
means_cm_df = cm_bins.groupby('Population_elastic_ID')[columns].mean().reset_index()

#### ---- Save the new mean df ---- ####
means_bp_df.to_csv(bp_output)
means_cm_df.to_csv(cm_output)




