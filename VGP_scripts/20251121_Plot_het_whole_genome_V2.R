#### Script to plot heterozygosity ####
## Date: 20251121 (November 21st, 2025)
## Author: Amanda Gardiner
## Version 2 (V1 is 20250325_Plot_het_whole_genome.R)
## GOAL: To plot heterozygosity calculated by 20250114_find_het_per_chr.sh 
##       I want to plot heterozygosity for each individual chromosome, save it, and then combine them for a whole genome graph
## NOTES: Want plots to follow those seen in Stanhope et al. 2023 Fig. 2
##        In Fig. 2A: Chromosomes on x axis across whole genome, and heterozygosity on y-axis
##        Het for each window in a chromosome will be plotted
##        In Fig. 2B: Histogram of windows across whole genome -- counting the number of windows with each value of heterozygosity
##        Modified from 20250123_Plot_het_per_chr.R
##        Made minor modifications to only save the whole genome plot, and not individual ones
## NOTES: Made V2 to simplify code, clean it up, and make sure it's plotting the tsv values correctly
####

#### ---- Import necessary libraries ---- ####
library(tidyverse)
library(ggplot2)
library(patchwork)

#### ---- Load in variables from shell ---- ####
args <- commandArgs()
dat <- read.table(args[6], header = FALSE)
colnames(dat) <- c(
  "Chr", "Start", "End", "Het", "Het_excl_ROH", 
  "Window_Size", "Window_Size_excl_ROH", 
  "Het_Per_KB", "Het_Per_KB_excl_ROH"
)
dat$Midpoint <- as.numeric(dat$Start) + (dat$Window_Size / 2)

chrom_length_file <- args[7]
output_file_name <- args[8]
num_aut_chr <- as.integer(args[9])
chroms <- unlist(strsplit(args[10:(num_aut_chr+10)], " "))

min_chrom_length <- 100000 ## Minimum length chromosomes have to be for calculating Het 

#### ---- Get chromosomes and filter for chroms > 100kb ---- ####
read_chromosome_file <- function(chrom_length_file) {
  dat <- read.table(chrom_length_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(dat) <- c("chr_name", "chr_length")
  return(dat)
}

chrom_df <- read_chromosome_file(chrom_length_file)
big_chroms <- chrom_df[chrom_df$chr_length > min_chrom_length & chrom_df$chr_name %in% chroms,"chr_name"]

## Order chromosomes
num_part <- as.numeric(gsub("[^0-9]", "", big_chroms))
chroms_ordered <- big_chroms[order(num_part)]

#### ---- Create base theme for each plot  ---- ####
base_theme <- theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(angle = 90, vjust = 0.5),
    plot.margin = margin(t = 0.5, r = 0.01, l = 0.25, b = 0.5)
  )

#### ---- Function to create plot by chromosome ---- ####
make_chr_plot <- function(idx, y_max) {
  chr_dat <- filter(dat, Chr == chroms_ordered[idx])
  if (nrow(chr_dat) == 0) return(NULL)
  
  fill_color <- ifelse(idx %% 2 == 0, "#707070", "#003C66")
  
  p <- ggplot(chr_dat, aes(x = Midpoint, y = Het_Per_KB_excl_ROH, group=1)) +
    geom_area(fill = fill_color) +
    base_theme +
    xlab(chr_dat$Chr) + 
    scale_y_continuous(limits = c(0, y_max), 
    trans="sqrt")
  
  # Add y-axis label only to the first chrom
  if (idx == 1) {
    p <- p + labs(y = "Heterozygosity per Kb")
  } else {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(t = 0.5, r = 0.01, l = 0.01, b = 0.5)
    )
  }
  
  return(p)
}
#### ---- END ---- ####

#### ---- PLOT, COMBINE, AND SAVE ---- ####

## Plot vesion 1, with all values shown and limit set around global max
global_max <- max(dat$Het_Per_KB_excl_ROH, na.rm = TRUE)
plot_list <- lapply(seq_along(chroms_ordered), function(i) make_chr_plot(i, global_max))

## Combine all plots and save
all_plots <- wrap_plots(plot_list, nrow = 1)
ggsave(output_file_name, 
  all_plots, 
  width = ifelse(length(plot_list)>=160, 100, length(plot_list)), 
  height = 10, 
  units = "in", 
  device = "png", 
  limitsize = FALSE)

