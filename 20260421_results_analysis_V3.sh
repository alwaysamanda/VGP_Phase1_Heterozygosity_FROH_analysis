#!/bin/bash
#### Script to create csv with all results and then run scripts to plot them ####
## Date: 20260421
## Version: 3 (V2 was 20251107_results_analysis_V2.sh)
## Author: Amanda Gardiner
## GOAL: Create a summary csv file of all the results and plot them
####

#### ---- Command for downloading the tetrapod trait database ---- ####
# ## Download tetrapod database
# wget https://zenodo.org/api/records/11303604/files/TetrapodTraits_1.0.0.csv/content
#### ---- END ---- ####


#### ---- Get today's date ---- ####
TODAY_DATE=$(date +%Y%m%d)


#### ---- Load in file names as variables ---- ####
CLADE_FILE=20251212_clades.txt
CLEANED_CLADE_FILE=20251212_clade_species_cleaned.txt
MAIN_FILE=../250430.VGP-Phase1/VGPPhase1-freeze-1.0.tsv
LINEAGE_FILE=20250930_VGP_Ordinal_List_Phase_1_Freeze.csv
PHYLO_FILE=250829.roadies_v1.1.4.nwk
MERGE_DATA="$TODAY_DATE"_FROH_Het_Trait_merged_data.csv
TETRAPOD_TRAITS=TetrapodTraits_1.0.0.csv
MISSING_LINEAGE_FILE=extended_lineage_missing_species.txt
CAPTIVITY_FILE=20260305_updated_captivity_flag.csv


#### ---- Get names of all cleaned  ---- ####
if [ -f "$CLEANED_CLADE_FILE" ]
    then
        echo "Cleaned clade file already made"
    else
        awk -F"/" '{print $3, $4}' $CLADE_FILE > $CLEANED_CLADE_FILE
        echo "Cleaned clade file finished"
fi


#### ---- Create database ---- ####
if [ -f "$MERGE_DATA" ]
    then
        echo "Database already created"
    else
        ## Run R script to get the data pulled and merged with trait data
        Rscript 20260217_create_database_V3.R $CLEANED_CLADE_FILE $TETRAPOD_TRAITS $MAIN_FILE $LINEAGE_FILE $MISSING_LINEAGE_FILE $CAPTIVITY_FILE
        echo "Finish create_database script"
        ## Add aligned genome length and ROH disparity between chromosomes to csv file
        python 20251107_Quantify_ROH_Patterns_V2.py $MERGE_DATA
        echo "Finish quantifying ROH patterns script"
fi