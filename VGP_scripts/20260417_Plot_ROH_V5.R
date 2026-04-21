#### Script to map ROH on chromosomes from ROH calculated by comparing assemblies ####
## Date: 20260417 (April 17th, 2026)
## Author: Amanda Gardiner
## Version 5 (Version 3 is 20250925_Plot_ROH_V4.R)
## GOAL: Plot ROH on each chromosome in the genome of a species
## NOTES: 
####

#### ---- Load in necessary packages ---- ####
library(tidyverse)
library(ggplot2)
library(grid)
library(svglite)
library(extrafont)
library(rphylopic)

#### ---- Load in data and variables ---- ####
args <- commandArgs()
Roh_file <- args[6]
chrom_length_file <- args[7]
aln_file <- args[8]
num_auto_chr <- as.numeric(args[9])
num_all_chr <- as.numeric(args[10])
genome_align_file <- args[11]
species_name <- args[12]
output_file <- args[13]

#### ---- Get Phylopic image for species ---- ####
get_image <- function(species_name){
  species_name_space <- str_to_sentence(str_replace(species_name, "_", " "))
  genus_name <- word(species_name_space, 1)
  
  for (name in c(species_name_space, genus_name)) {
    img <- tryCatch({
      uuid <- get_uuid(name = name, n = 1)
      if (name == genus_name) message("Using genus-level image for: ", species_name_space)
      get_phylopic(uuid = uuid)
    }, error = function(e) NULL)
    if (!is.null(img)) return(img)
  }
  message("Phylopic not found for: ", species_name_space, " — skipping image.")
  return(NULL)
}

#### ---- Read ROH data and add column to categorize based on ROH size ---- ####
read_roh_data <- function(roh_data_file){
    if (file.size(roh_data_file) == 0){
        dat <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
            c('Chr', 'Start', 'End', 'Length', 'ROH_Type'))
    }
    else{
        dat <- read.csv(roh_data_file, header=FALSE) ## Pull from the ROH_results.csv file
        colnames(dat) <- c('Chr', 'Start', 'End', 'Length')
        ROH_type <- c()
        for (i in 1:nrow(dat)) { 
            if (dat[i,4] >= 10000000) { 
                ROH_type[i] <- 'Ultralong' 
            } else if (dat[i,4] < 10000000 & dat[i,4] >= 1000000) {
                ROH_type[i] <- 'Long'
            } else if (dat[i,4] < 1000000 & dat[i,4] >= 500000) {
                ROH_type[i] <- 'Medium'
            } else {
                ROH_type[i] <- 'Short'
            } 
        }
        dat$ROH_Type <- ROH_type
        dat$ROH_Type <- factor(dat$ROH_Type, levels = c("Short", "Medium", "Long", "Ultralong"))
    }
    return(dat)
}

#### ---- Read in chromosomes and their lengths ---- ####
read_chrom_dat <- function(chrom_length_file, num_all_chr, num_auto_chr){
    Chrom <- read.table(chrom_length_file, header=FALSE)
    Chrom <- as.data.frame(Chrom[1:num_all_chr,])
    colnames(Chrom) <- c('Chrom_name', 'Chrom_end')
    return(Chrom)
}

#### ---- Read in alignment file ---- ####
read_aln_dat <- function(aln_file_name){
    dat <- read.table(aln_file_name, header=FALSE, sep="\t",stringsAsFactors = FALSE)
    colnames(dat) <- c("Chrom", "Start", "End")
    return(dat)
}

#### ---- Read in file with aligned genome proportion ---- ####
read_align_prop <- function(genome_align_file){
    lines <- readLines(genome_align_file)
    prop <- as.numeric(lines[4])
    return(prop)
}

#### ---- Plot results ---- ####
plot_ROH <- function(Roh_file, chrom_length_file, aln_file_name, num_all_chr, 
                     num_auto_chr, genome_align_file, species_name, file_name) {
  
  ## -- Data loading -- ##
  dat      <- read_roh_data(Roh_file)
  Chrom    <- read_chrom_dat(chrom_length_file, num_all_chr, num_auto_chr)
  aln_dat  <- read_aln_dat(aln_file_name)
  
  ## -- Label chromosomes and filter to autosomes only -- ##
  autosomal_names <- as.character(seq_len(num_auto_chr))
  sex_names       <- rep("Sex", num_all_chr - num_auto_chr)
  Chrom$Chrom_numbers <- c(autosomal_names, sex_names)
  
  ## -- Filter everything to autosomes only -- ##
  Chrom_auto  <- Chrom[Chrom$Chrom_numbers != "Sex", ]
  chrom_map   <- setNames(seq_len(num_auto_chr), autosomal_names)
  Chrom_auto$Chrom_pos <- chrom_map[Chrom_auto$Chrom_numbers]
  
  dat$Chrom_Label <- Chrom$Chrom_numbers[match(dat$Chr, Chrom$Chrom_name)]
  dat_auto        <- dat[dat$Chrom_Label != "Sex", ]
  dat_auto$Chrom_pos <- chrom_map[dat_auto$Chrom_Label]
  
  aln_dat$Chrom_Label <- Chrom$Chrom_numbers[match(aln_dat$Chrom, Chrom$Chrom_name)]
  aln_auto            <- aln_dat[aln_dat$Chrom_Label != "Sex", ]
  aln_auto$Chrom_pos  <- chrom_map[aln_auto$Chrom_Label]
  
  align_prop <- read_align_prop(genome_align_file)
  
  ## -- Get phylopic for species -- ##
  species_img <- get_image(species_name)

  ## --- Publication-ready theme --- ##
  pub_theme <- theme_classic(base_size = 11, base_family = "Myriad Web Pro") +
    theme(
      ## Axes
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.y     = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                       size = 8, colour = "black"),
      axis.text.y      = element_text(size = 8, colour = "black"),
      axis.title       = element_text(size = 10, face = "bold"),
      
      ## Legend
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.5, "cm"),
      legend.position  = "right",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      
      ## Panel
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      
      ## Annotation label
      plot.margin      = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt")
    )
  
  ## --- Build plot --- ##
  p <- ggplot() +
    
    ## Chromosome backbone
    geom_bar(
      data  = Chrom_auto,
      aes(x = Chrom_pos, y = Chrom_end),
      stat  = "identity",
      fill  = "grey85",
      colour = "grey70",
      width = 0.4
    ) +
    
    ## Alignment coverage (right offset)
    geom_segment(
      data = aln_auto,
      aes(x = Chrom_pos + 0.12, xend = Chrom_pos + 0.12,
          y = Start, yend = End),
      colour    = "#DB5461",
      linewidth = 0.8
    ) +
    
    ## ROH segments (left offset)
    geom_segment(
      data = dat_auto,
      aes(x = Chrom_pos - 0.12, xend = Chrom_pos - 0.12,
          y = Start, yend = End,
          colour = ROH_Type),
      linewidth = 1.2
    ) +
    
    ## Scales
    scale_x_continuous(
      breaks = chrom_map,
      labels = names(chrom_map),  
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_colour_manual(
      values = c(
        Ultralong = "#08519c",
        Long      = "#3182bd",
        Medium    = "#6baed6",
        Short     = "#bdd7e7"
      ),
      name  = "ROH class",
      guide = guide_legend(override.aes = list(linewidth = 2))
    ) +
    
    ## Genome alignment annotation
    annotate(
      "text",
      x     = max(chrom_map) + 0.4,
      y     = max(Chrom_auto$Chrom_end),
      label = paste0("Aligned: ", align_prop),
      hjust = 1, vjust = 1,
      size  = 3, colour = "grey30"
    ) +
    
    ## Labels
    labs(
      x = "",
      y = "Position (Mb)"
    ) +

    ## Phylopic
  {if (!is.null(species_img))
    add_phylopic(
      species_img,
      alpha = 1,
      x     = max(chrom_map) + 0.4,
      y     = max(Chrom_auto$Chrom_end) - (max(Chrom_auto$Chrom_end) / 4),
      ysize = max(Chrom_auto$Chrom_end) * 0.15)
  else
    NULL} +

    ## Theme
    pub_theme
  
  ## --- Save --- ##
    ggsave(
        filename = file_name,
        plot     = p,
        width    = 180,
        height   = 120,
        units    = "mm",
        device   = svglite
    )
}


plot_ROH(Roh_file, chrom_length_file, aln_file, num_all_chr, 
         num_auto_chr, genome_align_file, species_name, output_file)


