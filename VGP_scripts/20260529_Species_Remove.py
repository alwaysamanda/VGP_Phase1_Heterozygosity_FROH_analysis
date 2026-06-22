#### ==== Script for to list species to be removed from analysis for Gardiner et al. 2026 ==== ####
## Date: 20260529 (May 29th, 2026)
## Author: Amanda Gardiner
## Version: 1
## GOAL: Create a txt file listing all the species that have to be removed from the analysis for various reasons
## NOTES: 
#### ==== ==== ####


########## ========== LOAD IN NECESSARY PACKAGES ========== ##########
import sys
import os
import re
import numpy as np
import pandas as pd


########## ========== LOAD IN DATA AND VARIABLES ========== ##########
data                 = pd.read_csv(sys.argv[1])
output_species_list  = sys.argv[2]
output_filtered_data = sys.argv[3]


########## ========== CREATE SPECIES LIST ========== ##########
## -- Specific species that have come up with issues earlier in the alignment -- ##
species_list = [
    "Cetorhinus_maximus", ##Empty .1aln file
    "Girardinichthys_multiradiatus", ## Alternate turned out to be newer version of primary
    "Falco_punctatus", ## Very low alignment quality was messing with results
    "Micromys_minutus", ## Very low alignment quality was messing with results
    "Pristiophorus_japonicus", ## Empty .1aln file
    "Camelus_dromedarius" ## Something going on with the alignment -- getting a timeout error after 24hrs with variant calling
    "Heptrancihas_perlo" ## Something going on with the alignment -- It might be due to assembly error and VPG is checking
    "Sarcoramphus_papa" ## Something going on with alignment -- getting a timeout error after 24hrs with variant and alignment calling
]

## -- Species that have a very low alignment quality -- ##
data_aln_subset = data[data['genome_align_proportion'] <= .20] ## Remove all species which had less than 20% of their genome aligned
low_aln_list = list(data_aln_subset['species'])

## -- Species that had Heterozygosity calculated as 0 -- ##
data_no_het = data[data['heterozygosity'] == 0]
no_het_list = list(data_no_het['species'])

########## ========== COMPILE FINAL LIST AND SAVE IT ========== ##########
## -- Combine lists and remove duplicates -- ##
species_list.extend(x for x in low_aln_list if x not in species_list)
species_list.extend(x for x in no_het_list if x not in species_list)

## -- Save output file -- ##
with open(output_species_list, 'w') as output:
    output.write('\n'.join(species_list))


########## ========== REMOVE SPECIES FROM DF AND SAVE IT ========== ##########
## -- Filter data -- ##
data_fltr = data[~data['species'].isin(species_list)]

## -- Save the dataframe -- ##
data_fltr.to_csv(output_filtered_data, index=False)
