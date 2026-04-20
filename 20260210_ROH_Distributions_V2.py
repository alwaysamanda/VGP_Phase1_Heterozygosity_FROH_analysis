#### Script to plot SROH and NROH for different ROH categories ####
## Date: 20260210 (February 10th, 2026)
## Version: 2 (V1 is 20260119_ROH_Distributions.py)
## Author: Amanda Gardiner
## GOAL: Create plots which show the ROH distributions in a histogram/cumsum curve
##       Bins would be the different ROH size categories that we have
##       Will create a curve and save data for each species
## NOTES: Adding in function to calculate the aligned length of the ROHs
####

#### ---- Import necessary libraries ---- ####
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fpdf import FPDF
import re as re
import sys
from sys import stdin
import shutil
import io
import os

#### ----- Load in variables from the shell ---- ####
Aln_file = sys.argv[1]
num_aut_chr = int(sys.argv[2])
chromosome_list = sys.argv[3]
roh_data_file = sys.argv[4]
output_roh_gross_file = sys.argv[5]
output_roh_aln_file = sys.argv[6]
output_plot_gross_file = sys.argv[7]
output_plot_aln_file = sys.argv[8]
output_csv_file = sys.argv[9]
auto_chrom_names = sys.argv[10:(num_aut_chr+10)]

#### ---- FUNCTIONS ---- ####
## Function to read chromosome list
def read_chromosome_dat(chromosome_list_file, chromosome_names_list):
    Chrom = pd.read_csv(chromosome_list_file, header=None, sep='\t')
    Chrom.columns = ['Chrom_name', 'Chrom_end']
    if isinstance(chromosome_names_list, str):
            chromosome_names_list = [chromosome_names_list]
    Chrom_autosomal = Chrom[Chrom['Chrom_name'].isin(chromosome_names_list)]
    return Chrom_autosomal
## Finish function

## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[0] == chromosome]
    return dat
## Finish function

## Function to read alignment bed file
def aln_map_intervals(Alignment_list):
    starts = np.array([int(line[1]) for line in Alignment_list], dtype=np.int64) ## Start pos of each aln
    ends = np.array([int(line[2]) for line in Alignment_list], dtype=np.int64) ## End pos of each aln
    lengths = ends - starts
    
    # Sort intervals by start
    order = np.argsort(starts)
    starts, ends, lengths = starts[order], ends[order], lengths[order]
    
    # Get total aligned length
    total_length = np.sum(lengths)
    return total_length
## Finish function

## Define function to read the ROH data
def read_ROH_dat(roh_file_name):
    try:
        dat = pd.read_csv(roh_file_name, header=None)
        dat.columns = ['chrom', 'start', 'end', 'length']
        dat['length'] = pd.to_numeric(dat['length'], errors='coerce').astype('Int64')
    except pd.errors.EmptyDataError:
        print('Note: ROH csv was empty. Skipping.')
        dat = pd.DataFrame(columns=['chrom', 'start', 'end', 'length'])
    return dat
## Finish function

## Define function to categorize ROH into different categories
def categorize_ROH_vectorized(df, chrom_length, cum_lengths):
    bins = np.logspace(-3, 2, num=30, base=10.0) 
    ## For total length
    df["ROH_type"] = pd.cut((df["length"] / chrom_length) * 100, bins=bins, right=False)
    df["ROH_FROH_percent"] = (df["length"] / chrom_length) * 100
    ## Using ALN length for chrom and ROH
    df["ROH_type_aln"] = pd.cut((df["aln_length"] / cum_lengths) * 100, bins=bins, right=False)
    df["ROH_FROH_percent_aln"] = (df["aln_length"] / cum_lengths) * 100
    return df
## Finish function


## Function to find only the aligned regions within a ROH
def find_aln_ROH_regions(df, aln_file):
    aln_dat = pd.read_csv(aln_file, delimiter='\t')
    aln_dat.columns = ['chrom', 'start', 'end']
    results = []
    
    for _, row in df.iterrows():
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        
        ## Subset alignments for given chromosome
        subset_aln = aln_dat[aln_dat['chrom'] == chrom].copy()
        
        if len(subset_aln) == 0:
            results.append(0)
            continue
        
        ## Calculate overlap for each alignment segment
        overlap_start = np.maximum(subset_aln['start'].values, start)
        overlap_end = np.minimum(subset_aln['end'].values, end)
        overlaps = np.maximum(0, overlap_end - overlap_start)
        
        # Sum all overlaps
        total_overlap = overlaps.sum()
        results.append(total_overlap)
    
    return results
## Finish function


## Define function to calculate NROH/SROH/FROH for each chromosome and over the whole genome
def calc_ROH_stats(chromosome_list, Aln_file, roh_data_file, output_roh_gross_file, output_roh_aln_file, output_plot_gross_file, output_plot_aln_file, output_csv_file, auto_chrom_names):

    all_chroms = read_chromosome_dat(chromosome_list, auto_chrom_names)
    roh_dat = read_ROH_dat(roh_data_file)

    if roh_dat.empty:
        headers = pd.DataFrame(columns=["ROH_type", "NROH", "SROH", "FROH", "FROH_percent"])
        csv_headers = pd.DataFrame(columns=["Chrom", "Start", "End", "Length", "Aln_Length", "ROH_type", "FROH_percent", "ROH_type_aln", "FROH_percent_aln"])
        with open(output_roh_gross_file, "w") as file:
            headers.to_csv(file, sep="\t")
        with open(output_roh_aln_file, "w") as file:
            headers.to_csv(file, sep="\t")
        with open(output_csv_file, "w") as file:
            csv_headers.to_csv(file, sep=",")
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Arial", size=12)
        pdf.cell(200, 10, txt="No ROH were detected", ln=True, align='C')
        pdf.output(output_plot_gross_file)
        shutil.copy(output_plot_gross_file, output_plot_aln_file)
        return

    csv_results = []
    results_gross = []
    results_aln = []

    for _, row in all_chroms.iterrows():
        chrom = row["Chrom_name"]
        chrom_length = row["Chrom_end"]

        chrom_aln_list = read_aln_file(Aln_file, chrom)
        cum_lengths = aln_map_intervals(chrom_aln_list)

        chr_roh = roh_dat[roh_dat["chrom"] == chrom].copy()
        if chr_roh.empty:
            continue

        #### Calculate aligned length of ROH and % of chromosome ####
        chr_roh.insert(loc=4, column='aln_length', value='')
        chr_roh['aln_length'] = find_aln_ROH_regions(chr_roh, Aln_file)

        chr_roh = categorize_ROH_vectorized(chr_roh, chrom_length, cum_lengths)
        csv_results.append(chr_roh.reset_index())
        
        #### Calculate NROH and SROH ####
        grouped = chr_roh.groupby("ROH_type")["length"].agg(["count", "sum"]).rename(
            columns={"count": "NROH", "sum": "SROH"}
        )
        grouped["FROH"] = grouped["SROH"] / chrom_length
        grouped["FROH_percent"] = grouped["FROH"] * 100
        grouped["chrom"] = chrom
        results_gross.append(grouped.reset_index())
        
        grouped = chr_roh.groupby("ROH_type_aln")["aln_length"].agg(["count", "sum"]).rename(
            columns={"count": "NROH", "sum": "SROH"}
        )
        grouped["FROH"] = grouped["SROH"] / cum_lengths
        grouped["FROH_percent"] = grouped["FROH"] * 100
        grouped["chrom"] = chrom
        results_aln.append(grouped.reset_index())
        

    if not csv_results:
        csv_headers = pd.DataFrame(columns=["Chrom", "Start", "End", "Length", "Aln_Length", "ROH_type", "FROH_percent", "ROH_type_aln", "FROH_percent_aln"])
        with open(output_csv_file, "w") as file:
            csv_headers.to_csv(file, sep=",")

    if not results_gross:
        headers = pd.DataFrame(columns=["ROH_type", "NROH", "SROH", "FROH", "FROH_percent"])
        with open(output_roh_gross_file, "w") as file:
            headers.to_csv(file, sep="\t")
        with open(output_roh_aln_file, "w") as file:
            headers.to_csv(file, sep="\t")
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Arial", size=12)
        pdf.cell(200, 10, txt="No ROH were detected", ln=True, align='C')
        pdf.output(output_plot_gross_file)
        shutil.copy(output_plot_gross_file, output_plot_aln_file)
        return
    
    ## Concatenate all chromosomes
    results_df_gross = pd.concat(results_gross, ignore_index=True)
    results_df_aln = pd.concat(results_aln, ignore_index=True)
    csv_results_df = pd.concat(csv_results, ignore_index=True)

    ## Write results to file
    results_df_gross.to_csv(output_roh_gross_file, sep="\t", index=False)
    results_df_aln.to_csv(output_roh_aln_file, sep="\t", index=False)
    csv_results_df.to_csv(output_csv_file, sep=",", index=False)

    # Plot results
    grouped = results_df_gross.groupby("ROH_type")["FROH_percent"].sum()
    plt.bar(range(len(grouped)), grouped.values)
    plt.xticks(range(len(grouped)), [str(x) for x in grouped.index], rotation=45, ha='right')
    plt.ylabel("Sum of FROH_percent")
    plt.xlabel("ROH_type")
    plt.tight_layout() 
    plt.savefig(output_plot_gross_file)

    grouped = results_df_aln.groupby("ROH_type_aln")["FROH_percent"].sum()
    plt.bar(range(len(grouped)), grouped.values, alpha=0.75)
    plt.xticks(range(len(grouped)), [str(x) for x in grouped.index], rotation=45, ha='right')
    plt.ylabel("Sum of FROH_percent")
    plt.xlabel("ROH_type")
    plt.tight_layout() 
    plt.savefig(output_plot_aln_file)
## Finish function

#### ---- Run function ---- ####
calc_ROH_stats(chromosome_list, Aln_file, roh_data_file, output_roh_gross_file, output_roh_aln_file, output_plot_gross_file, output_plot_aln_file, output_csv_file, auto_chrom_names)
