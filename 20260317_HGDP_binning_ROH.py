#### Script to bin HGDP ROH ####
## Date: 20260317 (March 17th, 2026)
## Version: 1
## Author: Amanda Gardiner
## GOAL: Bin the HGDP roh based on size in cM and in bp
####

#### ---- Import necessary libraries ---- ####
import sys
from datetime import date
import numpy as np
import pandas as pd

#### ---- Load in data and variables ---- ####
args = sys.argv
date = date.today()
data = pd.read_csv(args[1])
output_bp_df = args[2]
output_cm_df = args[3]
log_breaks = np.logspace(-3, 2, 30)

#### ---- Create a df where they are binned in bp ---- ####
data['bp_log_bins'] = pd.cut(
    data['FROH_percent'],
    bins=log_breaks,
    include_lowest=True,
    right=False
)
data_raw = data.dropna(how='all')

bp_log_bin_summary = (
    data_raw
    .groupby(["IID", "Population_elastic_ID", "bp_log_bins"])
    .agg(FROH_percent=('FROH_percent', lambda x: x.sum(skipna=True)))
    .reset_index()
    .pivot(
        index=['IID', 'Population_elastic_ID'],
        columns='bp_log_bins',
        values='FROH_percent'
    )
    .fillna(0)
    .reset_index()
    .rename(columns={col: f'bp_log_bin_{col}' for col in data_raw['log_bins'].unique()}
            )
    .assign(log_bin_No_ROH=lambda x: (num_auto_chromosomes * 100) - x[
            [col for col in x.columns if col.startswith('bp_log_bin_')]
        ].sum(axis=1)
    )
)

#### ---- Save bp binned df ---- ####
bp_log_bin_summary.to_csv(output_bp_df)


#### ---- Create a df where they are binned in cM ---- ####
data['cm_log_bins'] = pd.cut(
    data['length_cm'],
    bins=log_breaks,
    include_lowest=True,
    right=False
)
data_raw = data.dropna(how='all')

cm_bin_summary = (
    data_raw
    .groupby(["IID", "Population_elastic_ID", "cm_log_bins"])
    .agg(length_cm=('length_cm', lambda x: x.sum(skipna=True)))
    .reset_index()
    .pivot(
        index=['IID', 'Population_elastic_ID'],
        columns='cm_log_bins',
        values='length_cm'
    )
    .fillna(0)
    .reset_index()
    .rename(columns={col: f'cm_log_bin_{col}' for col in data_raw['log_bins'].unique()}
            )
    .assign(log_bin_No_ROH=lambda x: (num_auto_chromosomes * 100) - x[
            [col for col in x.columns if col.startswith('cm_log_bin_')]
        ].sum(axis=1)
    )
)

#### ---- Save cM binned df ---- ####
cm_bin_summary.to_csv(output_cm_df)
