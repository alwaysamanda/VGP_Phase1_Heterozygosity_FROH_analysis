#### Script to create mask files for MSMC input to exclude ROH on any chromosome ####
## Date: 20250528 (May 28th, 2025)
## Author: Amanda Gardiner
## Version 1
## GOAL: Take ROH position txt files and use them to create negative mask files for MSMC input
## NOTES: Output file will be in bed format
##        ChromStart in column 2 -- (the first base on the chromosome is numbered 0 i.e. the number is zero-based)
##        ChromEnd in column 3 -- This position is non-inclusive, unlike chromStart (the first base on the chromosome is numbered 1 i.e. the number is one-based).

####

###############################################################################

## Import necessary libraries
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

## Load in variables from shell
clade = sys.argv[1]
spec_name = sys.argv[2]
chrom = sys.argv[3]
roh_dat = sys.argv[4]
output_file = sys.argv[5]

## FUNCTIONS ##
def READ_ROH_FILE(roh_file): 
    with open(roh_file, 'r') as file:
        header=next(file).strip().split(',')
        dat = [line.strip().split(',') for line in file if line.strip()]
    return pd.DataFrame(dat, columns=header) if dat else pd.DataFrame(columns=header)

def CREATE_NEGATIVE_BED_FILE(roh_file, output):
    dat = READ_ROH_FILE(roh_file)
    dat = dat.iloc[:, 0:3]
    dat.iloc[:, 1] = pd.to_numeric(dat.iloc[:, 1])
    dat.iloc[:, 2] = pd.to_numeric(dat.iloc[:, 2])
    dat.iloc[:, 2] = dat.iloc[:, 2] + 1 ## Adding 1 to each value since the end value is not included in mask
    dat.to_csv(output, sep='\t', index=False, header=False)

CREATE_NEGATIVE_BED_FILE(roh_dat, output_file)





