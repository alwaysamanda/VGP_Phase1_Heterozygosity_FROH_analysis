# Inbreeding and Heterozygosity in VGP Phase 1 species
Author: Amanda Gardiner and Richard Durbin

First: April 20th, 2026

## Overview
This repository provides the scripts and codes used in the paper by Gardiner et al. (In Prep) for 'Genetic diversity, runs of homozygosity, and demographic history across vertebrates is associated with extinction risk.'

The software here is used to take FastGA-aligned assemblies and calculate the length and positions of Runs of Homozygosity (ROH), heterozygosity outside of ROH across the genome, and $N_e$ (effective population size) estimated using the sequential markovian coalescent.  This was done on the VGP Phase 1 assemblies, with all primaries coming from the VGP while some secondaries were supplemented. 

This work was performed using resources provided by the Cambridge Service for Data Driven Discovery (CSD3) operated by the University of Cambridge Research Computing Service (www.csd3.cam.ac.uk), provided by Dell EMC and Intel using Tier-2 funding from the Engineering and Physical Sciences Research Council (capital grant EP/T022159/1), and DiRAC funding from the Science and Technology Facilities Council (www.dirac.ac.uk).

## Table of Contents
- [Input_data](#Input_data)
- [VGP_analysis](#VGP_analysis)
- [HGDP_analysis](#HGDP_analysis)
- [Citation](#Citation)
- [Contact](#Contact)
- [References](#References)


## Input_data


## VGP_analysis

Snakemake was used to automate the analysis for all species.  Conda environments were used for running the snakefile and also for all downstream scripts. All conda packages used are in the file Conda_packages_list.txt.

All code that was used to run analysis for each individual species in available in the direcotry VGP_scripts. Each species had a config file created, which are all stored in config_files subdirectory, which is organized by clade and then shortname for each species (first letter of clade followed by first three letters in genus and species each, for example Homo sapiens becomes mHomSap). The profiles subdirectory contains the config.yaml file which was used to specify the base paraemters for running snakemake on an HPC.

The snakefile for running snakemake is Snakefile_filtering. All custom scripts written for the pipeline are included in the directory. The scripts mask_file.py and indel_mask.py were written by Jaskaran Singh Gill and originally found [here](https://github.com/vicegill/Cambridge-Primate-Analysis).

Additionally required are the script paftools.js from minimap2, which is available [here](https://github.com/lh3/minimap2/tree/master/misc); the script block_bootstrap_mhs.py from cobraa by Trevor Cousins, available [here](https://github.com/trevorcousins/cobraa); and the scripts generate_multihetsep.py and bamCaller.py from msmc-tools, available [here](https://github.com/stschiff/msmc-tools).

Once analysis is finished for all species, the script 20260421_results_analysis_V3.sh is used to created a csv with data for all species.


## HGDP_analysis

Snakemake was used to automate the analysis for the HGDP data. The snakefile used (snakefile_hgdp) and all dependent scripts are found in the HGDP_scripts directory, along with downstream analyses in the HGDP_scripts/HGDP_analyses. All conda packages used are in the file Conda_packages_list.txt.
HGDP vcf files and mask file were downloaded from https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/



## Citation
To be added

## Contact
Amanda Gardiner: ag2427@cam.ac.uk

## References
Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. “Sustainable Data Analysis with Snakemake.” 10:33. Preprint, F1000Research, April 19, 2021. https://doi.org/10.12688/f1000research.29032.2.


