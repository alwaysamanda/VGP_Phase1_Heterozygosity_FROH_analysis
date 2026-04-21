#### Script to get calculate ROH using RZooRoH ####
## Date: 20260318 (March 18th, 2026)
## Version: 1
## Author: Amanda Gardiner
## GOAL: Calculate ROH in centimorgans for HGDP data using RZooRoH package
####

#### ---- Import necessary libraries ---- ####
library(RZooRoH)
library(tidyverse)

#### ---- Load in data and variables ---- ####
args <- commandArgs(trailingOnly = TRUE)
date <- Sys.Date()
gen <- read.table(args[1], header=FALSE, stringsAsFactors=FALSE)
cat("=== INPUT GEN CLASSES ===\n")
print(sapply(gen, class))
cat("=== INPUT GEN HEAD ===\n")
print(head(gen, 3))
gen$V2 <- gsub(":", "_", gen$V2)
cm_file <- read.table(args[2], header = FALSE)
colnames(cm_file) <- c('snp', 'chrom', 'pos')
output_gen_file <- args[3]
output_results <- args[4]
output_segments_csv <- args[5]

## Check if gen has duplicated SNP label columns and remove duplicate if necessary
if (identical(gen$V2, gen$V3)) {
  gen <- subset(gen, select = -V3)
}

stopifnot(nrow(gen) == nrow(cm_file))

#### ---- Convert SNP positions to centimorgans ---- ####
gen <- left_join(gen, cm_file %>% select(snp, pos), by = c("V2" = "snp")) %>%
    mutate(V4 = pos) %>%
    select(-pos)

con <- gzfile(output_gen_file, "w")
write.table(gen, con, quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)

check <- read.table(gzfile(output_gen_file), header=FALSE, nrows=5)
cat("=== DIMENSIONS ===\n")
print(dim(check))
cat("=== COLUMN CLASSES ===\n")
print(sapply(check, class))
cat("=== HEAD ===\n")
print(head(check))
cat("=== UNIQUE VALUES IN GENO COLS (first 3 individuals) ===\n")
print(lapply(check[, 5:7], unique))

#### ---- Run RZooROH to get ROH positions on chromosome ---- ####
data <- zoodata(
        genofile = output_gen_file,
        zformat = "gt",
        min_maf = 0.01
    )

model <- zoomodel(
        K = 4,
        nT = snakemake@threads, 
        base = 0.5,     
        err = 0.005     
        )

results <- zoorun(
        zoomodel = model,
        zooin = data,
        nT = 4,
        localhbd = FALSE      
        )

roh_segments <- results@hbdseg

#### ---- Save the output ---- ####
saveRDS(results, output_results)
write.csv(roh_segments, output_segments_csv, row.names = FALSE)