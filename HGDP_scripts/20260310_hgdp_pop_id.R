#### Script to assign population IDs to the samples from HGDP ####
## Date: 202560310 (March 10th, 2026)
## Version: 1
## Author: Amanda Gardiner
## GOAL: Assign population IDs to samples from HGDP
####

#### ---- Load in necessary packages ---- ####
library(tidyverse)

#### ---- Load in data and variables ---- ####
args <- commandArgs()
date <- Sys.Date()
data <- read_csv(args[6])
pop_info_file <- read_tsv(args[7])
output_file_name <- args[8]

#### ---- Merge files based on sample ID numbers ---- ####
pop_dat <- subset(pop_info_file, select = c("Sample_name", "Population_elastic_ID"))
names(pop_dat)[names(pop_dat) == "Sample_name"] <- "IID"

data <- data %>%
  left_join(pop_dat, by = "IID")

#### ---- Save output file ---- ####
write.csv(data, output_file_name)