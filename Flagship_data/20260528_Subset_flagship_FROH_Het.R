#### ==== Script for plotting Het/ROH figures for Formenti et al. 2026 ==== ####
## Date: 20260528 (May 28th, 2026)
## Author: Amanda Gardiner
## GOAL: Create jitterplots of heterozygosity outside of ROH regions and ROH for VGP Phase 1 species
## NOTES: 
####

########## ========== Load in necessary packages ========== ##########
library(ggplot2)
library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(see)
library(ggpubr)
library(ggpmisc)
library(systemfonts)
library(svglite)

########## ========== Load in data and variables ========== ##########
args <- commandArgs()
date <- Sys.Date()
data <- read_csv(args[6])
results_directory <- args[7]


########## ========== Filtering criteria ========== ##########
## -- Filter out DD/NE species for plotting and factor the remaining species -- ##
data <- data %>% filter(IUCN != "DD/NE", !is.na(IUCN))
iucn_levels <- c("CR", "EN", "VU", "NT", "LC")
data$IUCN <- factor(data$IUCN, levels=iucn_levels)

## -- Filter out invertebrates -- ##
data <- data %>% filter(data$Extended_lineage != "Cyclostomes" & data$Extended_lineage != "Other Deuterostomes" )

## -- Filter out for wild only species -- ##
data_wild <- data %>% filter(zoo_flag == "no")

## -- Remove species which had erroneous results for Het calculations -- ##
data_wild <- data_wild %>% filter(heterozygosity > 0)


########## ========== Set themes for all plots ========== ##########
My_Theme = theme(
  text = element_text(size = 30, family = "Myriad Web Pro"),
  axis.text = element_text(size = 20, hjust = 0.5, family = "Arial"),
  axis.title.x = element_text(size = 20, hjust = 0.5, family = "Arial"),
  axis.title.y = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(), 
  legend.position = "inside",
  legend.position.inside = c(1, 0.8),
  legend.justification = c(1, 0.5),
  axis.line = element_line(colour = 'black', linewidth = 1), 
  axis.ticks = element_blank(),
  axis.ticks.length=unit(.25, "cm"),
  plot.margin = margin(5, 15, 5, 5)
  )

#### ---- Set IUCN color palette ---- ####
iucn_colors <- c(
    "LC" = "#60C659", 
    "NT" = "#CCE226", "VU" = "#F9E814", 
    "EN" = "#FC7F3F", "CR" = "#D81E05" 
  )

########## ========== MAIN FIGURE ========== ##########
#### ---- Heterozygosity for wild species ---- ####
het_wild <- ggplot(data_wild, aes(x=IUCN, y=heterozygosity, color=IUCN, fill=IUCN)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
  ) +
  geom_jitter(
    position = position_jitter(width = 0.15, seed = 43),
    shape = 21,
    color = "black",
    aes(fill = IUCN),
    stroke = 1,
    size = 8,
    alpha = 0.8
  ) +
  xlab('') +
  ylab('Heterozygosity per kb outside ROH ≥100kb') +
  scale_fill_manual(
    values = iucn_colors,
    drop = FALSE,
    breaks = c("CR", "EN", "VU", "NT", "LC")
  ) +
  guides(color = guide_legend(title = "Group")) +
  scale_y_continuous(
    trans = scales::trans_new(
        name      = "log1p_base10",
        transform = function(x) log10(pmax(x, 0) + 1),
        inverse   = function(x) 10^x - 1,
        breaks    = scales::extended_breaks(),
        domain    = c(0, Inf)
    ),
    breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 30),
    labels = c("0.01", "0.05", "0.1", "0.5", "1", "2", "5", "10", "20", "30"),
    expand = expansion(mult = c(0.01, 0.03)),
    limits = c(-0.05, NA),
    oob = scales::squish
    ) + 
  My_Theme +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  coord_flip(clip = "off")
#### ---- END ---- ####

#### ---- FROH for ROH ≥1% of their chromosome, in wild data ---- ####
p_froh_longplus_wild <- ggplot(
  data_wild, 
  aes(
    x=IUCN, 
    y=`1%+_aln_FROH_percent`, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
  ) +
  geom_jitter(
    position = position_jitter(width = 0.15, seed = 43),
    color = "black",
    aes(fill = IUCN), 
    shape = 21,
    stroke = 1,
    size=8, 
    alpha = 0.8
  ) + 
  scale_fill_manual(
    values = iucn_colors, 
    drop = FALSE
  ) +  
  guides(color = guide_legend(title = "Group")) + 
  My_Theme + 
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome",
    trans = scales::trans_new(
        name    = "log1p_base10",
        transform = function(x) log10(pmax(x, 0) + 1),
        inverse   = function(x) 10^x - 1,
        breaks    = scales::extended_breaks(),
        domain    = c(0, Inf)
    ),
    breaks = c(0, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
    labels = c("0", "1", "2", "5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"),
    expand = expansion(mult = c(0.01, 0.03)),
    limits = c(-0.05, NA), 
    oob = scales::squish
    ) + 
  theme(
    legend.position = "none"
  ) + 
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  coord_flip(clip = "off") 
#### ---- END ---- ####


#### ---- Assemble Figure ---- ####
Fig_one <- (
    p_froh_longplus_wild / 
    het_wild
)

#### ---- Save Figure ---- ####
Fig_one_name <- paste0(results_directory, "Main_Figures/", date, "_Subset_Main_Figure_1.svg") 
ggsave(
  Fig_one_name,
  plot = Fig_one,
  width = 10, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
########## ========== END FIGURE ========== ##########