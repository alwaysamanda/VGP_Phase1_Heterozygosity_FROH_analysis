# Inbreeding and Heterozygosity in VGP Phase 1 species
Author: Amanda Gardiner and Richard Durbin

First: April 20th, 2026

## Overview
This repository provides the scripts and codes used in the paper by Gardiner et al. (In Prep) for 'Genetic diversity, runs of homozygosity, and demographic history across vertebrates is associated with extinction risk.'

The software here is used to take FastGA-aligned assemblies and calculate the length and positions of Runs of Homozygosity (ROH), heterozygosity outside of ROH across the genome, and $N_e$ (effective population size) estimated using the sequential markovian coalescent.  This was done on the VGP Phase 1 assemblies, with all primaries coming from the VGP while some secondaries were supplemented. 


## Table of Contents
- [Input data and genome alignment](#Input data)
- [VGP analysis](#VGP analysis)
- [HGDP analysis](#HGDP analysis)


## Input data


## VGP analysis

Snakemake was used to automate the analysis for all species. The snakefile used (Snakefile_filtering) and all dependent scripts are found in the VGP_scripts directory. Conda environments were used for running the snakefile and also for downstream plotting. All conda packages used are in the file Conda_packages_list.txt.

### ROH analysis
We custom wrote a software to calculate ROH when comparing two assemblies. We calculate ROH based on the equation:
$$
S_{i} = S_{i-1} + (X_{i} - X_{i-1}) - P
$$
where $S_{i}$ is the length `score' for ROH at any given variant, $S_{i-1}$ is the score at the previous variant, $(X_{i} - X_{i-1})$ represents the distance in aligned bases between the current variant and previous variant, and $P$ is the penalty threshold that must be passed to be considered an ROH. 


## HGDP analysis

Snakemake was used to automate the analysis for the HGDP data. The snakefile used (snakefile_hgdp) and all dependent scripts are found in the HGDP_scripts directory, along with downstream analyses in the HGDP_scripts/HGDP_analyses. All conda packages used are in the file Conda_packages_list.txt.
HGDP vcf files and mask file were downloaded from https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/



## Citation
To be added

## Contact
Amanda Gardiner: ag2427@cam.ac.uk
