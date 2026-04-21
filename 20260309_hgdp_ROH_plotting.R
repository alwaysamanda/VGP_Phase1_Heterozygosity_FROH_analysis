#### Script to plot ROH patterns from HGDP data ####
## Date: 202560309 (March 9th, 2026)
## Version: 1
## Author: Amanda Gardiner
## GOAL: Visualize distribution of ROH in HGDP data
####

#### ---- Load in necessary packages ---- ####
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpmisc)
library(ggrepel)
library(ggridges)
library(gplots)
library("RColorBrewer")
library(viridis)

#### ---- Load in data and variables ---- ####
args <- commandArgs()
date <- Sys.Date()
data <- read_csv(args[6])
chrom_file <- read_delim(args[7], col_names=FALSE, delim=' ')
colnames(chrom_file) <- c("Chrom", "length")
num_auto_chromosomes <- 22


#### ---- Set themes and colors for plotting ---- ####
eps <- 1e-3
lilcol <- viridis(100)
bigcol <- viridis(200)

pop_colors <- c(
  "AdygeiHGDP" = "#284600", 
  "BalochiHGDP" = "#d8ffdf", 
  "BantuKenyaHGDP" = "#f2c2ff", 
  "BantuSouthAfricaHGDP" = "#f2c2ff", 
  "BasqueHGDP" = "#284600", 
  "BedouinHGDP" = "#4e2100", 
  "BergamoItalianHGDP" = "#284600", 
  "BiakaHGDP" = "#f2c2ff", 
  "BougainvilleHGDP" = "#001991", 
  "BrahuiHGDP" = "#d8ffdf", 
  "BurushoHGDP" = "#d8ffdf", 
  "CambodianHGDP" = "#b35200", 
  "ColombianHGDP" = "#ff9081", 
  "DaiHGDP" = "#b35200", 
  "DaurHGDP" = "#b35200", 
  "DruzeHGDP" = "#4e2100", 
  "FrenchHGDP" = "#284600", 
  "HanHGDP" = "#b35200", 
  "HazaraHGDP" = "#d8ffdf", 
  "HezhenHGDP" = "#b35200", 
  "JapaneseHGDP" = "#b35200", 
  "KalashHGDP" = "#d8ffdf", 
  "KaritianaHGDP" = "#ff9081", 
  "LahuHGDP" = "#b35200", 
  "MakraniHGDP" = "#d8ffdf", 
  "MandenkaHGDP" = "#f2c2ff", 
  "MayaHGDP" = "#ff9081", 
  "MbutiHGDP" = "#f2c2ff", 
  "MiaoHGDP" = "#b35200", 
  "MongolianHGDP" = "#b35200", 
  "MozabiteHGDP" = "#4e2100", 
  "NaxiHGDP" = "#b35200", 
  "NorthernHanHGDP" = "#b35200", 
  "OrcadianHGDP" = "#284600", 
  "OroqenHGDP" = "#b35200", 
  "PalestinianHGDP" = "#4e2100", 
  "PapuanHighlandsHGDP" = "#001991", 
  "PapuanSepikHGDP" = "#001991", 
  "PathanHGDP" = "#d8ffdf", 
  "PimaHGDP" = "#ff9081", 
  "RussianHGDP" = "#284600", 
  "SanHGDP" = "#f2c2ff", 
  "SardinianHGDP" = "#284600", 
  "SheHGDP" = "#b35200", 
  "SindhiHGDP" = "#d8ffdf", 
  "SuruiHGDP" = "#ff9081", 
  "TuHGDP" = "#b35200", 
  "TujiaHGDP" = "#b35200", 
  "TuscanHGDP" = "#284600", 
  "UygurHGDP" = "#d8ffdf", 
  "XiboHGDP" = "#b35200", 
  "YakutHGDP" = "#b35200", 
  "YiHGDP" = "#b35200", 
  "YorubaHGDP" = "#f2c2ff"
)

#### ---- Functions ---- ####
extract_bin_start <- function(bin_name) {
  bin_name <- sub("log_bin_", "", bin_name)
  start_val <- as.numeric(sub("\\[([0-9.e+-]+),.*", "\\1", bin_name))
  return(start_val)
}

#### ---- Bin and aggregate data ---- ####
min_length <- min(data$length[data$length > 0], na.rm = TRUE)
max_length <- max(data$length, na.rm = TRUE)

data_long <- data %>%
    mutate(
        length = as.numeric(length),
        FROH_percent = as.numeric(FROH_percent),
        IID = IID,
        Pop_ID = Population_elastic_ID,
        log_bins_percent = cut(
          FROH_percent, 
          breaks = 10^seq(-3, log10(max(FROH_percent, na.rm = TRUE)), length.out = 30),
          include.lowest = TRUE, 
          right = FALSE),
        log_bins_length = cut(
          length,
          breaks = 10^seq(log10(min_length), log10(max_length), length.out = 30),
          include.lowest = TRUE,
          right = FALSE), 
        linear_bins_percent = cut(
          FROH_percent, 
          breaks = seq(0, 100, by = 2),
          include.lowest = TRUE, 
          right = FALSE),
        sqrt_bins_percent = cut(
          FROH_percent, 
          breaks = seq(0, sqrt(100), length.out = 30)^2,
          include.lowest = TRUE, 
          right = FALSE)
      )

log_bin_summary <- data_long %>%
  group_by(IID, Population_elastic_ID, log_bins_percent) %>%
  summarise(FROH_percent = sum(FROH_percent, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = log_bins_percent, 
              values_from = FROH_percent,
              names_prefix = "log_bin_",
              values_fill = 0)

heatmap_log_dat <- log_bin_summary %>%
  select(-Population_elastic_ID) %>%
  column_to_rownames(var = "IID")

col_names <- colnames(heatmap_log_dat)
roh_bin_cols <- grep("^log_bin_", col_names, value = TRUE)
other_cols <- setdiff(col_names, roh_bin_cols)
bin_starts <- sapply(roh_bin_cols, extract_bin_start)
sorted_roh_bins <- roh_bin_cols[order(bin_starts)]
heatmap_log_dat <- heatmap_log_dat[, c(sorted_roh_bins, other_cols)]

## L1 normalize rows
row_sums <- rowSums(heatmap_log_dat) # heatmap_log_dat[, sorted_roh_bins]
heatmap_log_dat <- heatmap_log_dat / ifelse(row_sums == 0, 1, row_sums)


## Assign population to each individual
indiv_names <- unique(log_bin_summary$IID)
group_colors <- pop_colors[log_bin_summary$Population_elastic_ID]
names(group_colors) <- log_bin_summary$IID
group_colors <- group_colors[!duplicated(names(group_colors))]
row_color_bar <- group_colors[indiv_names]

#### ---- Do a heatmap of binned data (1 row per individual) ---- ####
png(paste0("Figures/", date, "_HGDP_Log_bin_heatmap.png"), width = 3000, height = 6000,
    units = "px", res=300)
heatmap.2(as.matrix(heatmap_log_dat), 
          scale = "none", 
          col = lilcol, 
          trace = "none", 
          density.info = "none",  
          dendrogram = "row", 
          hclustfun = function(d) hclust(d, method = "average"),
          Colv = FALSE, 
          margins = c(12, 15), 
          cexCol = 0.4, 
          cexRow = 0.4, 
          RowSideColors = row_color_bar, 
          key.par = list(mar = c(1, 0.5, 3, 1)), 
          keysize = 0.5
          )
dev.off()

## Do heatmaps with non-ROH bins included
total_roh_per_species <- data %>%
  filter(!is.na(FROH_percent)) %>%  # Exclude the "No_ROH" placeholder rows
  group_by(IID, Population_elastic_ID) %>%
  summarise(
    total_ROH_percent = sum(FROH_percent, na.rm = TRUE),
    .groups = "drop"
  )

non_roh_summary <- total_roh_per_species %>%
  mutate(
    total_genome_percent = num_auto_chromosomes * 100, 
    non_ROH_percent = total_genome_percent - total_ROH_percent
  ) %>%
  select(IID, Population_elastic_ID, non_ROH_percent)

data_read_log <- log_bin_summary %>%
  left_join(non_roh_summary, by = c("IID", "Population_elastic_ID"))

heatmap_log_dat <- data_read_log %>%
  select(-Population_elastic_ID) %>% 
  column_to_rownames(var = "IID")

col_names <- colnames(heatmap_log_dat)
roh_bin_cols <- grep("^log_bin_", col_names, value = TRUE)
other_cols <- setdiff(col_names, roh_bin_cols)
bin_starts <- sapply(roh_bin_cols, extract_bin_start)
sorted_roh_bins <- roh_bin_cols[order(bin_starts)]
heatmap_log_dat <- heatmap_log_dat[, c(sorted_roh_bins, other_cols)]

## L1 normalize rows
row_sums <- rowSums(heatmap_log_dat) # heatmap_log_dat[, sorted_roh_bins]
heatmap_log_dat <- heatmap_log_dat / ifelse(row_sums == 0, 1, row_sums)


## Assign population to each individual
indiv_names <- unique(log_bin_summary$IID)
group_colors <- pop_colors[log_bin_summary$Population_elastic_ID]
names(group_colors) <- log_bin_summary$IID
group_colors <- group_colors[!duplicated(names(group_colors))]
row_color_bar <- group_colors[indiv_names]

log_breaks <- eps * (1/eps) ^ seq(0, 1, length.out = 201)

png(paste0("Figures/", date, "_HGDP_Log_all_heatmap.png"), width = 3000, height = 6000,
    units = "px", res=300)
heatmap.2(as.matrix(heatmap_log_dat), 
          scale = "none",
          breaks = log_breaks, 
          col = bigcol, 
          trace = "none", 
          density.info = "none",  
          dendrogram = "row", 
          hclustfun = function(d) hclust(d, method = "average"), 
          Colv = FALSE, 
          margins = c(12, 15), 
          cexCol = 0.4, 
          cexRow = 0.4, 
          RowSideColors = row_color_bar, 
          key.par = list(mar = c(1, 1, 1, 1)),  
          keysize = 0.5,  
          key.xlab = "Log scaled prop.", 
        )
dev.off()


#### ---- Do a heatmap but with ROH in bases and not normalized FROH percent ---- ####
log_bin_summary_bases <- data_long %>%
  group_by(IID, Population_elastic_ID, log_bins_length) %>%
  summarise(length = sum(length, na.rm = TRUE), .groups = "drop") %>%
  select(-Population_elastic_ID) %>%                              # drop before complete()
  complete(IID, log_bins_length, fill = list(length = 0)) %>%
  pivot_wider(names_from = log_bins_length, 
              values_from = length,
              names_prefix = "log_bin_",
              values_fill = 0)

heatmap_log_dat <- log_bin_summary_bases %>%
  select(-matches("log_bin_NA")) %>%                              # drop NA bin if present
  column_to_rownames(var = "IID")

col_names <- colnames(heatmap_log_dat)
roh_bin_cols <- grep("^log_bin_", col_names, value = TRUE)
other_cols <- setdiff(col_names, roh_bin_cols)
bin_starts <- sapply(roh_bin_cols, extract_bin_start)
sorted_roh_bins <- roh_bin_cols[order(bin_starts)]
heatmap_log_dat <- heatmap_log_dat[, c(sorted_roh_bins, other_cols)]

write.csv(heatmap_log_dat, "HGDP_ROH_Bins.csv", row.names=TRUE)

## Assign population to each individual
indiv_names <- unique(log_bin_summary$IID)
group_colors <- pop_colors[log_bin_summary$Population_elastic_ID]
names(group_colors) <- log_bin_summary$IID
group_colors <- group_colors[!duplicated(names(group_colors))]
row_color_bar <- group_colors[indiv_names]

#### ---- Do a heatmap of binned data (1 row per individual) ---- ####
png(paste0("Figures/", date, "_HGDP_SROH_Log_bin_heatmap.png"), width = 3000, height = 6000,
    units = "px", res=300)
heatmap.2(as.matrix(heatmap_log_dat), 
          scale = "none", 
          col = lilcol, 
          trace = "none", 
          density.info = "none",  
          dendrogram = "row", 
          hclustfun = function(d) hclust(d, method = "average"),
          Colv = FALSE, 
          margins = c(12, 15), 
          cexCol = 0.4, 
          cexRow = 0.4, 
          RowSideColors = row_color_bar, 
          key.par = list(mar = c(1, 0.5, 3, 1)), 
          keysize = 0.5
          )
dev.off()

png(paste0("Figures/", date, "_HGDP_SROH_correlationdist_Log_bin_heatmap.png"), width = 3000, height = 6000,
    units = "px", res=300)
heatmap.2(as.matrix(heatmap_log_dat), 
          scale = "none", 
          col = lilcol, 
          trace = "none", 
          density.info = "none",  
          dendrogram = "row", 
          distfun = function(x) as.dist(1 - cor(t(x))), 
          hclustfun = function(d) hclust(d, method = "average"),
          Colv = FALSE, 
          margins = c(12, 15), 
          cexCol = 0.4, 
          cexRow = 0.4, 
          RowSideColors = row_color_bar, 
          key.par = list(mar = c(1, 0.5, 3, 1)), 
          keysize = 0.5
          )
dev.off()
#### ---- END ---- ####

# #### ---- Do a ROH map of human chromosomes but make it a heatmap by presence of ROH ---- ####
# n_bins <- 200000 

# ## Specify the lengths of the chromosomes -- taken from GRCH38
# chr_lengths <- data.frame(
#   chrom_names = c(
#     1, 2, 3, 
#     4, 5, 6, 
#     7, 8, 9, 
#     10, 11, 12, 
#     13, 14, 15, 
#     16, 17, 18, 
#     19, 20, 21, 
#     22), 
#   lengths = c(
#     248956422, 242193529, 198295559,
#     190214555, 181538259, 170805979,
#     159345973, 145138636, 138394717,
#     133797422, 135086622, 133275309,
#     114364328, 107043718, 101991189,
#     90338345,  83257441,  80373285,
#     58617616,  64444167,  46709983,
#     50818468)
# )

# ## Calculate the density of the ROH across chromosomes
# chr_length_vec <- setNames(chr_lengths$lengths, chr_lengths$chrom_names)
# binned <- data %>%
#   group_by(CHR) %>%
#   group_modify(~ {
#     chr_name <- as.character(.y$CHR)
#     chr_len  <- chr_length_vec[[chr_name]]     
#     breaks   <- seq(0, chr_len, length.out = n_bins + 1)
#     bin_mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
#     bin_size <- chr_len / n_bins

#     .x %>%
#       rowwise() %>%
#       reframe(
#         y_bin = which(breaks[-length(breaks)] < POS2 & breaks[-1] > POS1),
#       ) %>%
#       mutate(
#         y_mid    = bin_mids[y_bin],
#         bin_size = bin_size
#       ) %>%
#       group_by(y_bin, y_mid, bin_size) %>%
#       summarise(count = n(), .groups = "drop")
#   }) %>%
#   ungroup()

# binned_complete <- binned %>%
#   group_by(CHR) %>%
#   complete(y_bin = seq(min(y_bin), max(y_bin), by = 1),
#            fill = list(count = 0)) %>%
#   ungroup()

# p <- ggplot(binned_complete, aes(x = factor(CHR), y = y_bin, fill = count)) +
#   geom_tile(height = 1, width = 10) +
#   scale_fill_viridis_c(option = "magma", name = "Segment\nCount") +
#   scale_y_continuous(expand = c(0, 0)) +
#   facet_wrap(~ CHR, scales = "free_y", nrow = 1, strip.position = "bottom") +
#   labs(y = "Position") +
#   theme_minimal() +
#   theme(
#     axis.title.x  = element_blank(),
#     axis.text.x   = element_blank(),
#     axis.ticks.x  = element_blank(),
#     axis.text.y   = element_blank(),
#     axis.ticks.y  = element_blank(),
#     panel.grid    = element_blank(),
#     panel.spacing = unit(0.1, "lines")
#   )


# chrom_heatmap_figure_name = paste0("Figures/", date, "_HGDP_ROH_chrom_heatmap.png")
# png(file = chrom_heatmap_figure_name, width = 3000, height = 2000, res = 300)
# print(p)
# dev.off()