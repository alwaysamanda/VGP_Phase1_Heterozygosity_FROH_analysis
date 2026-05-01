#### Script for plotting figures for paper on VGP ROH and Het data####
## Date: 20260420 (April 20th, 2026)
## Author: Amanda Gardiner
## Version: 3
## GOAL: Create main and supplementary figures for publication
## NOTES: 
####

#### ---- Load in necessary packages ---- ####
library(ggplot2)
library(tidyverse)
library(data.table)
library(patchwork)
library(magick)
library(grid)
library(see)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(ggsignif)
library(svglite)
library(systemfonts)
library(dunn.test)
library(MASS)
library(colorspace)
library(vegan)
library(hexbin)
library(plotly)
library(geometry)
library(ggridges)
library(rstatix)


#### ---- Load in data and variables ---- ####
args <- commandArgs()
date <- Sys.Date()
data <- read_csv(args[6])
HGDP_dat <- read_csv(args[7])
HGDP_autosome_length <- as.numeric(2900000000)
results_directory <- args[8]


#### ---- Filtering criteria for data ---- ####
## -- Filter out DD/NE species for plotting and factor the remaining species -- ##
data <- data %>% filter(IUCN != "DD/NE", !is.na(IUCN))
data$IUCN <- factor(data$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))


## -- Filter out invertebrates -- ##
data <- data %>% filter(data$Extended_lineage != "Cyclostomes" & data$Extended_lineage != "Other Deuterostomes" )


## -- Normalize NROH by the number of chromosomes -- ##
data <- data %>% mutate(NROH_normalized = all_aln_NROH / num_auto_chromosomes)


## -- Factor the extended lineages -- ##
data$Extended_lineage <- factor(
    data$Extended_lineage, 
    levels=c(
        "Mammals", 
        "Birds", 
        "Crocodilians", 
        "Turtles", 
        "Lepidosauria", 
        "Amphibians", 
        "Lobe-finned fishes", 
        "Ray-finned fishes", 
        "Cartilaginous fishes"
    ))


## -- Filter out for wild only species -- ##
data_wild <- data %>% filter(Combined_captivity_flag_primary == "no" | Combined_captivity_flag_primary == "no likely")

## -- Set colors for HGDP populations -- ##
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

#### ---- Set theme for all plots ---- ####
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


#### ---- Set habitat color palettes ---- ####
Marine_flag_colors <- c(
    "yes" = "#1A7A8A", 
    "no" = "#C8C5BC"
  )

microhabitat_colors <- c(
  "Ter" = "#C7A96A", 
  "Fos" = "#A0522D", 
  "Arb" = "#4CA66B", 
  "Aqu" = "#2196A6",
  "Aer" = "#6BAED6"
)




########## ========== MAIN FIGURE 1 ========== ##########
#### ---- Heterozygosity for wild species ---- ####
het_wild <- ggplot(data_wild, aes(x=IUCN, y=heterozygosity, color=IUCN, fill=IUCN)) +
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1) + 
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN), 
    stroke = 1,
    size=8, 
    alpha = 0.8
  ) + 
  xlab('') + 
  ylab('Heterozygosity per kb outside ROH ≥100kb') + 
  scale_fill_manual(values = iucn_colors, 
      drop = FALSE, 
      breaks = c("CR", "EN", "VU", "NT", "LC") 
      ) +  
  guides(color = guide_legend(title = "Group")) + 
  scale_y_continuous(
    trans = "sqrt", 
    breaks = seq(0, 30, by=2),
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  My_Theme + 
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  coord_flip() 
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
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    },
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1
  ) + 
  geom_jitter(
    position = position_jitter(0.15),
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
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  theme(
    legend.position = "none"
  ) + 
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  coord_flip() 
#### ---- END ---- ####



#### ---- Assemble Figure 1 ---- ####
Fig_one <- (
    het_wild /
    p_froh_longplus_wild
)


## Save figure
Fig_one_name <- paste0(results_directory, "Main_Figures/", date, "_Main_Figure_1_mean.svg") 
ggsave(
  Fig_one_name,
  plot = Fig_one,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
#### ---- END ---- ####




########## ========== MAIN FIGURE 2 ========== ##########
#### ---- Subset and prepare wild data ---- ####
dat_clean_wild <- data_wild %>%
  filter(!is.na(IUCN)) %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized))


## -- Set axis limits and plot reference line where 1 NROH is one long ROH -- ##
x_max <- max(dat_clean_wild$all_aln_FROH_percent, na.rm = TRUE)
y_max <- max(dat_clean_wild$NROH_normalized, na.rm = TRUE)
 
x_buffer <- x_max * 0.05
y_buffer <- y_max * 0.05

x_min_plot <- 0
x_max_plot <- x_max + x_buffer
y_min_plot <- 0
y_max_plot <- y_max + y_buffer

y_max_ref <- x_max_plot


## -- Get convex hulls -- ##
roh_hulls <- dat_clean_wild %>%
  group_by(IUCN) %>%
  filter(n() >= 3) %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  slice(chull(all_aln_FROH_percent, NROH_normalized))

## -- Get density-weighted centroids -- ##
weighted_centroids <- dat_clean_wild  %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  group_by(IUCN) %>%
  group_modify(~ {
    kde <- MASS::kde2d(
      .x$all_aln_FROH_percent,
      .x$NROH_normalized,
      n = 100
    )
    x_idx <- findInterval(.x$all_aln_FROH_percent, kde$x)
    y_idx <- findInterval(.x$NROH_normalized, kde$y)
    x_idx <- pmax(1, pmin(x_idx, length(kde$x) - 1))
    y_idx <- pmax(1, pmin(y_idx, length(kde$y) - 1))
    weights <- kde$z[cbind(x_idx, y_idx)]
    weights <- weights / sum(weights)
    data.frame(
      cx = sum(.x$all_aln_FROH_percent * weights),
      cy = sum(.x$NROH_normalized * weights)
    )
  }) %>%
  ungroup()

#### ---- Plot for wild data ---- ####
p_overlay <- ggplot(dat_clean_wild, aes(x = all_aln_FROH_percent, y = NROH_normalized,
                                    color = IUCN, fill = IUCN)) +
  annotate("segment",
           x = 0, xend = x_max_plot,
           y = 0, yend = y_max_ref,
           linetype = "dashed", color = "grey30",
           linewidth = 0.7, alpha = 0.8) +
  geom_polygon(
    data = roh_hulls,
    aes(fill = IUCN, color = IUCN),
    alpha = 0.10,
    linewidth = 0.7
  ) +
  geom_point(
    data = weighted_centroids,
    aes(x = cx, y = cy, fill = IUCN),
    size = 6,
    shape = 21,
    color = "black",
    stroke = 1.2,
    alpha = 0.75,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = iucn_colors, guide = "none") +
  scale_fill_manual(
    values = iucn_colors,
    name = "IUCN status",
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        alpha = 0.75,
        stroke = 1.2,
        color = "black"
      )
    )
  ) +
  scale_x_sqrt(
    limits = c(x_min_plot, x_max_plot),
    expand = c(0, 0),
    breaks = c(0, 5, 10, 20, 30, 50, 75, 100),
    labels = c("0", "5", "10", "20", "30", "50", "75", "100")
  ) +
  scale_y_sqrt(
    limits = c(y_min_plot, y_max_plot),
    expand = c(0, 0),
    breaks = c(0, 100, 500, 1000, 2000, 4000, 6000, 8000),
    labels = c("0", "100", "500", "1000", "2000", "4000", "6000", "8000")
  ) +
  labs(
    x = expression(bold(F[ROH] ~ "(%)" ~ (sqrt ~ scale))),
    y = expression(bold(N[ROH]))
  ) +
  My_Theme + 
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey60",
                                     linewidth = 0.4),
    legend.key.size = unit(0.55, "cm"),
    legend.text = element_text(size = 9, color = "grey10", family = "Arial"),
    legend.title = element_text(size = 10, face = "bold", color = "grey10",
                                family = "Arial"),
    plot.margin = margin(0, 0, 5, 5)
  )

p_x <- ggplot(dat_clean_wild, aes(x = all_aln_FROH_percent, y = IUCN,
                              fill = IUCN, color = IUCN)) +
  geom_jitter(aes(color = IUCN),
              height = 0.18, size = 2.8, alpha = 0.55, shape = 21,
              stroke = 0.6) +
  geom_density_ridges(
    alpha = 0.6,
    linewidth = 0.7,
    scale = 1.3,
    quantile_lines = TRUE,
    quantiles = 2,
    color = "black",
    panel_scaling = TRUE
  ) +
  scale_fill_manual(values = iucn_colors, guide = "none") +
  scale_color_manual(values = iucn_colors, guide = "none") +
  scale_x_sqrt(
    limits = c(x_min_plot, x_max_plot),
    expand = c(0, 0),
    breaks = c(0, 5, 10, 20, 30, 50, 75, 100),
    labels = c("0", "5", "10", "20", "30", "50", "75", "100")
  ) +
  labs(x = NULL, y = "IUCN status") +
  My_Theme + 
  theme(
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "grey20", family = "Arial"),
    plot.margin = margin(5, 0, 0, 5)
  )

p_y <- ggplot(dat_clean_wild, aes(x = NROH_normalized, y = IUCN,
                              fill = IUCN, color = IUCN)) +
  geom_jitter(aes(color = IUCN),
              height = 0.18, size = 2.8, alpha = 0.55, shape = 21,
              stroke = 0.6) +
  geom_density_ridges(
    alpha = 0.6,
    linewidth = 0.7,
    scale = 1.3,
    quantile_lines = TRUE,
    quantiles = 2,
    color = "black",
    panel_scaling = TRUE
  ) +
  scale_fill_manual(values = iucn_colors, guide = "none") +
  scale_color_manual(values = iucn_colors, guide = "none") +
  scale_x_sqrt(
    limits = c(y_min_plot, y_max_plot),
    expand = c(0, 0),
    breaks = c(0, 100, 500, 1000, 2000, 4000, 6000, 8000),
    labels = c("0", "100", "500", "1000", "2000", "4000", "6000", "8000")
  ) +
  coord_flip() +
  labs(x = expression(bold(N[ROH])), y = NULL) +
  My_Theme + 
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "grey20", family = "Arial"),
    plot.margin = margin(0, 5, 5, 0)
  )

layout <- "
AAAAAAB
CCCCCCD
"

p_combined <- p_x + plot_spacer() + p_overlay + p_y +
  plot_layout(
    design = layout,
    widths = c(5, 1),
    heights = c(1, 2.5)
  ) +
  plot_annotation(
    title = "ROH distribution across IUCN conservation categories",
    subtitle = "Convex hulls with density-weighted centroids and marginal distributions",
    theme = theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 14, face = "bold", color = "grey10",
                                hjust = 0, family = "Arial",
                                margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0,
                                   face = "italic", family = "Arial",
                                   margin = margin(b = 8)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 12, 12, 12)
    )
  )


Fig_two_name <- paste0(results_directory, "Main_Figures/", date, "_Main_Figure_2.svg")
ggsave(
  Fig_two_name, 
  p_combined,
  width = 18, 
  height = 14, 
  device = svglite)

#### ---- Run statistical tests for wild data ---- ####
stats_wild_output_file <- paste0(results_directory, "Main_Figures/", date, "_Fig_2_stats.txt")

## -- Prepare scaled data -- ##
dat_ord <- dat_clean_wild %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  mutate(
    froh_sqrt = sqrt(all_aln_FROH_percent),
    nroh_sqrt = sqrt(NROH_normalized),
    froh_scaled = scale(froh_sqrt),
    nroh_scaled = scale(nroh_sqrt)
  )

dist_mat <- dist(cbind(dat_ord$froh_scaled, dat_ord$nroh_scaled))

## -- PERMANOVA -- ##
set.seed(42)
permanova_result <- adonis2(
  dist_mat ~ IUCN,
  data = dat_ord,
  permutations = 999,
  by = "margin"
)

## -- Betadisper + permutest -- ##
betadisp_result <- betadisper(dist_mat, dat_ord$IUCN)
betadisp_test <- permutest(betadisp_result, permutations = 999)
betadisp_tukey <- TukeyHSD(betadisp_result)

## -- Pairwise PERMANOVA -- ##
iucn_levels <- levels(dat_ord$IUCN)
pairs <- combn(iucn_levels, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(pair) {
  idx <- dat_ord$IUCN %in% pair
  sub_dist <- as.dist(as.matrix(dist_mat)[idx, idx])
  sub_data <- dat_ord[idx, ]
  set.seed(42)
  res <- adonis2(sub_dist ~ IUCN, data = sub_data, permutations = 999)
  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    R2 = round(res$R2[1], 3),
    F_stat = round(res$F[1], 3),
    p_value = res$`Pr(>F)`[1]
  )
})

pairwise_df <- do.call(rbind, pairwise_results) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "bonferroni"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

## -- Convex hull areas -- ##
hull_areas <- dat_clean_wild %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  group_by(IUCN) %>%
  filter(n() >= 3) %>%
  summarise(
    n = n(),
    hull_area_raw = {
      h <- chull(all_aln_FROH_percent, NROH_normalized)
      pts <- cbind(all_aln_FROH_percent[h], NROH_normalized[h])
      x <- pts[, 1]
      y <- pts[, 2]
      abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
    },
    hull_area_sqrt = {
      h <- chull(sqrt(all_aln_FROH_percent), sqrt(NROH_normalized))
      pts <- cbind(sqrt(all_aln_FROH_percent)[h], sqrt(NROH_normalized)[h])
      x <- pts[, 1]
      y <- pts[, 2]
      abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
    },
    .groups = "drop"
  ) %>%
  arrange(desc(hull_area_raw))

## -- Write all results to file -- ##
sink(stats_wild_output_file)
cat("=== PERMANOVA (centroid differences) ===\n")
print(permanova_result)
cat("\n=== Betadisper (dispersion/hull size differences) ===\n")
print(betadisp_test)
cat("\n=== Pairwise Tukey HSD on dispersion ===\n")
print(betadisp_tukey)
cat("\n=== Pairwise PERMANOVA (Bonferroni corrected) ===\n")
print(pairwise_df)
cat("\n=== Convex hull areas per IUCN group ===\n")
print(hull_areas)
sink()

cat("Stats written to:", stats_wild_output_file, "\n")
#### ---- END FOR WILD DATA ---- #####



########## ========== MAIN FIGURE 3 ========== ##########
#### ---- Set variables for MSMC ---- ####
mu <- 1.25e-8

#### ---- Read MSMC data ---- ####
data_read <- map2_dfr(
  data_wild$MSMC_file_name,
  data_wild$species,
  ~{
    if (!file.exists(.x)) {
      message("Skipping missing file: ", .x)
      return(NULL)
    }
    spec_dat <- readr::read_table(.x, col_names = TRUE) %>%
  mutate(right_time_boundary = as.numeric(right_time_boundary))
    tibble(
      Generations_end   = spec_dat$left_time_boundary / mu,
      Generations_start = spec_dat$right_time_boundary / mu, 
      mid_gen           = (spec_dat$left_time_boundary + spec_dat$right_time_boundary) / (2 * mu),
      Ne                = (1 / spec_dat$lambda) / (2 * mu),
      Species = .y
    )
  }
) %>%
  left_join(dplyr::select(data, species, IUCN),
            by = c("Species" = "species"))
data_read$IUCN <- factor(data_read$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))


## -- Build common log time bins -- ##
time_grid <- exp(seq(
  log(min(data_read$Generations_start[data_read$Generations_start > 0])),
  log(max(data_read$Generations_end)),
  length.out = 200
))

## -- Interpolate all species onto the bins
data_interp <- data_read %>%
  group_by(Species, IUCN) %>%
  group_modify(~{
    tibble(
      Generations = time_grid,
      Ne = approx(
        x    = .x$mid_gen,
        y    = .x$Ne,
        xout = time_grid,
        rule = 2
      )$y
    )
  }) %>%
  ungroup()


#### ---- Summarize data and save file ---- ####
data_summary <- data_interp %>%
  group_by(IUCN, Generations) %>%
  summarise(
    median_Ne = median(Ne, na.rm = TRUE),
    lq_Ne     = quantile(Ne, 0.25, na.rm = TRUE),
    uq_Ne     = quantile(Ne, 0.75, na.rm = TRUE),
    .groups   = "drop"
  )

write.csv(data_summary, paste0(results_directory, date, "_wild_MSMC_grouped_stats.csv"))


#### ---- Plot MSMC grouped by IUCN status ---- ####
wild_iucn <- ggplot(
  data_summary,
  aes(
    x     = Generations, 
    y     = median_Ne, 
    color = IUCN, 
    fill  = IUCN)) +
  geom_ribbon(
    aes(ymin = lq_Ne, ymax = uq_Ne), 
    alpha = 0.3, 
    colour = NA) + 
  geom_line(linewidth = 1) +
  scale_color_manual(values = iucn_colors, drop = FALSE) +  
  scale_fill_manual(values = iucn_colors, drop = FALSE) + 
  scale_x_log10(labels = scales::label_scientific()) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(
    x = "Generations",
    y = "Effective Population Size",
    color = "IUCN Status",
    fill  = "IUCN Status"
  ) +
  My_Theme + 
  annotation_logticks(sides = "bl")
#### ---- END ---- ####


#### ---- Jitterplot of contemoprary NE across IUCN status ---- ####
modern_dat <- data_read[which(data_read$Generations_end== 0), ]
modern_dat$IUCN <- factor(modern_dat$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))

p_modern_ne <- ggplot(
  modern_dat, 
  aes(
    x=IUCN, 
    y=Ne, 
    fill=IUCN
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1
  ) + 
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN), 
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
  theme(
    legend.position = "none"
  ) + 
  scale_y_continuous(
    trans = "log10", 
    expand = expansion(mult = c(0.01, 0.03))
  ) +  
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  annotation_logticks(sides = "b") + 
  coord_flip()
#### ---- END ---- ####

#### ---- Assemble Figure 3 ---- ####
Fig_three <- (
    wild_iucn /
    p_modern_ne
)


## Save figure
Fig_three_name <- paste0(results_directory, "Main_Figures/", date, "_Main_Figure_3_Mean.svg") 
ggsave(
  Fig_three_name,
  plot = Fig_three,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
#### ---- END ---- ####




########## ========== MAIN FIGURE 4 ========== ##########
#### ---- Get habitat information for species ---- ####
data_wild <- rename(data_wild, c(longplus_froh=`1%+_aln_FROH_percent`))

wild_gen_hab_dat <- data_wild[c(
  "species", 
  "IUCN", 
  "Combined_captivity_flag_primary", 
  "habitats", 
  "longplus_froh", 
  "all_aln_FROH_percent", 
  "heterozygosity"
  )]

wild_gen_hab_dat <- wild_gen_hab_dat %>%
  mutate(
    major_habitat_count = case_when(
      is.na(habitats) ~ NA_integer_,
      TRUE ~ habitats %>%
        str_split(";\\s*") %>%
        lapply(\(x) str_extract(x, "^[^_]+")) %>%
        sapply(\(x) length(unique(x))) %>%
        as.integer()
    ),
    all_habitat_count = case_when(
      is.na(habitats) ~ NA_integer_,
      TRUE ~ lengths(str_split(habitats, ";\\s*"))
    )
  )

data_subset_wild <- data_wild[c("species", "IUCN", "longplus_froh", "all_aln_FROH_percent", "heterozygosity", "Fos", "Aer", "Aqu", "Ter", "Arb")]
data_subset_long_wild <- data_subset_wild %>%
  pivot_longer(
    cols = 6:10, 
    names_to = "group", 
    values_to = "presence"
  )
data_subset_long_wild <- data_subset_long_wild %>% filter(presence == 1)

## Set order for microhabitats in plots
data_subset_long_wild$group <- factor(data_subset_long_wild$group, levels=c("Aer", "Arb", "Fos", "Ter", "Aqu"))

#### ---- Plot species which occur in marine environments versus those that don't with FROH for ROH ≥1% of their chromosome in length ---- ####
df_major_hab <- wild_gen_hab_dat %>%
  mutate(habitats = trimws(habitats),
         Marine_flag = if_else(str_detect(habitats, "(^|; )(9|10|11|12)_"), "yes", "no"))

df_major_hab <- df_major_hab %>% filter(!is.na(Marine_flag))


## -- Test for statistically significant differences in FROH values first -- ##
stat_test <- df_major_hab %>%
  rstatix::dunn_test(longplus_froh ~ Marine_flag, p.adjust.method = "BH")

sig_pairs <- stat_test %>%
  filter(p.adj < 0.05) %>%
  mutate(annotation = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*"
  ))

## -- Plot -- ##
marine_flag_wild_froh_long <- ggplot(
  data = df_major_hab, 
  aes(
    x = Marine_flag, 
    y = longplus_froh, 
    fill=Marine_flag)) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.5, 
    lwd = 1.1) + 
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 1,
    size=8, 
    alpha = 0.8
  ) + 
  scale_fill_manual(values=Marine_flag_colors
  ) +
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  coord_flip() + 
  My_Theme

## -- Test for statistically significant differences in FROH by marine flag -- ##
kruskal.test(longplus_froh ~ Marine_flag, data = df_major_hab)


#### ---- Plot microhabitat against FROH for ROH ≥1% of their chromosome in length ---- ####
## -- Test for statistically significant differences in FROH values first -- ##
stat_test <- data_subset_long_wild %>%
  rstatix::dunn_test(longplus_froh ~ group, p.adjust.method = "BH")

sig_pairs <- stat_test %>%
  filter(p.adj < 0.05) %>%
  mutate(annotation = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*"
  ))

comparisons_list <- Map(c, sig_pairs$group1, sig_pairs$group2)
annotations_list <- sig_pairs$annotation

## -- Plot -- ##
microhabitat_plot_wild_frohlong <- ggplot(
  data_subset_long_wild, 
  aes(
    x = group, 
    y = longplus_froh
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.6,
    color = "black", 
    lwd = 1.1, 
    aes(fill=group)
    ) + 
  geom_jitter(
    alpha = 0.5, 
    position = position_jitter(
      width = 0.1, 
      height = 0), 
    size = 8, 
    shape = 21, 
    color = "black", 
    aes(fill=group)
  ) +
  geom_signif(
    comparisons = comparisons_list,
    annotations = annotations_list,
    step_increase = 0.04,    
    tip_length = 0.01,
    textsize = 5
  ) +
  coord_flip(clip = "off") + 
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by = 10), 
    labels = as.character(seq(0, 100, by = 10)),
    expand = expansion(mult = c(0.01, 0.01)), 
    limits = c(0, 105)
  ) + 
  scale_fill_manual(values = microhabitat_colors) + 
  My_Theme + 
  theme(plot.margin = margin(t = 5, r = 80, b = 5, l = 5, unit = "pt"))


#### ---- Microhabitat for Heterozygosity ---- ####
## Test for statistically significant differences in het values
stat_test <- data_subset_long_wild %>%
  rstatix::dunn_test(heterozygosity ~ group, p.adjust.method = "BH")

sig_pairs <- stat_test %>%
  filter(p.adj < 0.05) %>%
  mutate(annotation = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*"
  ))

comparisons_list <- Map(c, sig_pairs$group1, sig_pairs$group2)
annotations_list <- sig_pairs$annotation

## Plot
microhabitat_het <- ggplot(
  data_subset_long_wild, 
  aes(
    x = group, 
    y = heterozygosity
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.6,
    color = "black", 
    lwd = 1.1, 
    aes(fill=group)
    ) + 
  geom_jitter(
    alpha = 0.5, 
    position = position_jitter(
      width = 0.1, 
      height = 0), 
    size = 8, 
    shape = 21, 
    color = "black", 
    aes(fill=group)
  ) +
  geom_signif(
    comparisons = comparisons_list,
    annotations = annotations_list,
    step_increase = 0.04,    
    tip_length = 0.01,
    textsize = 5
  ) +
  coord_flip(clip = "off") + 
  scale_y_continuous(
    trans = "sqrt", 
    breaks = seq(0, 30, by=2),
    expand = expansion(mult = c(0.01, 0.01)), 
    limits = c(0, 47)
  ) +  
  scale_fill_manual(
    values=microhabitat_colors
  ) + 
  My_Theme + 
  theme(plot.margin = margin(t = 5, r = 80, b = 5, l = 5, unit = "pt"))

## -- Test for statistically significant differences in het between microhabitats -- ##
kruskal.test(heterozygosity ~ group, data = data_subset_long_wild)

dunn.test::dunn.test(
  data_subset_long_wild$heterozygosity,
  data_subset_long_wild$group,
  method = "bh"
)

#### ---- END ---- ####


#### ---- Assemble Figure 4 ---- ####
Fig_four <- (
    marine_flag_wild_froh_long /
    microhabitat_plot_wild_frohlong /
    microhabitat_het
) +
  plot_layout(
    guides = "collect",
    axes = "collect"
  ) &
  theme(
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 80, b = 5, l = 5, unit = "pt")
  )


## Save figure
Fig_four_name <- paste0(results_directory, "Main_Figures/", date, "_Main_Figure_4.svg") 
ggsave(
  Fig_four_name,
  plot = Fig_four,
  width = 30, 
  height = 15,
  units = "in",
  dpi = 600,
  device=svglite 
)
#### ---- END ---- ####










########## ========== SUPPLEMENTARY FIGURES ========== ##########
#### ---- Heterozygosity for all species ---- ####
het_alldat <- ggplot(
  data, 
  aes(x=IUCN, y=heterozygosity, color=IUCN, fill=IUCN)) +
  geom_boxplot(
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1) + 
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN), 
    stroke = 1,
    size=8, 
    alpha = 0.8
  ) + 
  xlab('') + 
  ylab('Heterozygosity per kb outside ROH ≥100kb') + 
  scale_fill_manual(values = c(
      "LC" = "#60C659", 
      "NT" = "#CCE226", "VU" = "#F9E814", 
      "EN" = "#FC7F3F", "CR" = "#D81E05"
      ), 
      drop = FALSE, 
      breaks = c("CR", "EN", "VU", "NT", "LC") 
      ) +  
  guides(color = guide_legend(title = "Group")) + 
  scale_y_continuous(
    trans = "sqrt", 
    breaks = seq(0, 30, by=2),
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  My_Theme + 
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  coord_flip() 

## Save figure
Het_alldat_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_AllData_Het_Fig.svg") 
ggsave(
  Het_alldat_figname,
  plot = het_alldat,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300,
  device=svglite 
)
#### ---- END ---- ####



#### ---- FROH for ROH ≥100kb, in all data ---- ####
p_frohall_alldat <- ggplot(
  data_wild, 
  aes(
    x=IUCN, 
    y=`all_aln_FROH_percent`, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1
  ) + 
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN), 
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
    name = "Inbreeding Coefficient (%) for ROH ≥100kb of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  theme(
    legend.position = "none"
  ) + 
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  coord_flip() 
#### ---- END ---- ####

## Save figure
froh_all_alldat_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_Wild_FROH_all_Fig.svg") 
ggsave(
  froh_all_alldat_figname,
  plot = p_frohall_alldat,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300,
  device=svglite 
)
#### ---- END ---- ####




#### ---- FROH for ROH ≥1% of their chromosome, in all data ---- ####
p_froh_longplus_alldat <- ggplot(
  data, 
  aes(
    x=IUCN, 
    y=`1%+_aln_FROH_percent`, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1
  ) + 
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN), 
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
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  theme(
    legend.position = "none"
  ) + 
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  coord_flip() 


## Save figure
froh_long_alldat_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_AllData_FROH_long_Fig.svg") 
ggsave(
  froh_long_alldat_figname,
  plot = p_froh_longplus_alldat,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300,
  device=svglite 
)
#### ---- END ---- ####



#### ---- FROH for ROH ≥100kb, in all data ---- ####
p_frohall_alldat <- ggplot(
  data, 
  aes(
    x=IUCN, 
    y=`all_aln_FROH_percent`, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1
  ) + 
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN), 
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
    name = "Inbreeding Coefficient (%) for ROH ≥100kb of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  theme(
    legend.position = "none"
  ) + 
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  coord_flip() 


## Save figure
froh_all_alldat_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_AllData_FROH_all_Fig.svg") 
ggsave(
  froh_all_alldat_figname,
  plot = p_frohall_alldat,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300,
  device=svglite 
)
#### ---- END ---- ####




#### ---- HGDP NROH/SROH plots ---- ####
roh_summary <- HGDP_dat %>%
  group_by(IID, Population_elastic_ID) %>%
  summarise(
    ROH_count = n(),
    ROH_median = median(length),
    ROH_max_value = max(length), 
    ROH_length_sum = sum(length),
    .groups        = "drop"
  )

roh_summary_clean <- roh_summary %>%
  filter(!is.na(Population_elastic_ID))

## -- Compute weighted centroids per population -- ##
## Weight = number of individuals in population (n)
pop_centroids <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  summarise(
    n               = n(),
    centroid_sroh   = weighted.mean(ROH_length_sum, w = rep(1, n())),  # equal weights = mean
    centroid_nroh   = weighted.mean(ROH_count,      w = rep(1, n())),
    sd_sroh         = sd(ROH_length_sum),
    sd_nroh         = sd(ROH_count),
    se_sroh         = sd(ROH_length_sum) / sqrt(n()),
    se_nroh         = sd(ROH_count)      / sqrt(n()),
    .groups         = "drop"
  )

## Compute hulls
roh_hulls <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  filter(n() >= 3) %>%
  slice(chull(ROH_length_sum, ROH_count))

## Order populations by centroid position (SROH then NROH)
pop_order <- pop_centroids %>%
  arrange(centroid_sroh, centroid_nroh) %>%
  pull(Population_elastic_ID)

roh_summary_clean <- roh_summary_clean %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

roh_hulls <- roh_hulls %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

pop_centroids <- pop_centroids %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

## -- Get limits of plots -- ##
x_min <- min(roh_summary_clean$ROH_length_sum)
x_max <- max(roh_summary_clean$ROH_length_sum)
y_min <- min(roh_summary_clean$ROH_count)
y_max <- max(roh_summary_clean$ROH_count)

## -- Plot reference lines for Long ROH = 1Mb -- ##
y2_max <- x_max / 1000000

#### ---- PLOT ---- ####
p <- ggplot(roh_summary_clean, aes(x = ROH_length_sum, y = ROH_count,
                                    color = Population_elastic_ID,
                                    fill  = Population_elastic_ID)) +
  annotate("segment",
           x = 0, xend = x_max,
           y = 0, yend = y2_max,
           linetype = "longdash", color = "grey40", linewidth = 0.4) +
  geom_point(size = 1.5, alpha = 0.4) +
  geom_polygon(data = roh_hulls, alpha = 0.2, linewidth = 0.3) +
  ## -- Weighted centroid: larger point + crosshairs (SE) -- ##
  geom_errorbar(
    data = pop_centroids,
    aes(x = centroid_sroh,
        ymin = centroid_nroh - se_nroh,
        ymax = centroid_nroh + se_nroh),
    width = 0, linewidth = 0.5, inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbarh(
    data = pop_centroids,
    aes(y    = centroid_nroh,
        xmin = centroid_sroh - se_sroh,
        xmax = centroid_sroh + se_sroh),
    height = 0, linewidth = 0.5, inherit.aes = FALSE,
    color = "black"
  ) +
  geom_point(
    data = pop_centroids,
    aes(x = centroid_sroh, y = centroid_nroh, fill = Population_elastic_ID),
    shape = 21, size = 3, color = "black", stroke = 0.6, inherit.aes = FALSE
  ) +
  scale_color_manual(values = pop_colors) +
  scale_fill_manual(values  = pop_colors) +
  facet_wrap(~ Population_elastic_ID) +
  labs(
    x = "SROH (sum of ROH length)",
    y = "NROH (number of ROH)"
  ) +
  theme_classic() +
  theme(
    strip.text    = element_text(size = 6),
    axis.text     = element_text(size = 5),
    axis.title    = element_text(size = 8),
    legend.position = "none"
  )

HGDP_fig_name <- paste0(results_directory, "Supplementary_Figures/", date, "_HGDP_pop_hulls.svg")
ggsave(
  HGDP_fig_name,
  p,
  width  = 20,
  height = 16,
  dpi    = 600,
  device = svglite
)
#### ---- END ---- ####




#### ---- VGP NROH/SROH plots for all data ---- ####
#### ---- Subset and prepare all data ---- ####
dat_clean <- data %>%
  filter(!is.na(IUCN)) %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized))
 
## Set axis limits and plot reference land
x_max <- max(dat_clean$all_aln_FROH_percent, na.rm = TRUE)
y_max <- max(dat_clean$NROH_normalized, na.rm = TRUE)
 
x_buffer <- x_max * 0.05
y_buffer <- y_max * 0.05

x_min_plot <- 0
x_max_plot <- x_max + x_buffer
y_min_plot <- 0
y_max_plot <- y_max + y_buffer

y_max_ref = x_max_plot

## -- Get convex hulls -- ##
roh_hulls <- dat_clean %>%
  group_by(IUCN) %>%
  filter(n() >= 3) %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  slice(chull(all_aln_FROH_percent, NROH_normalized))

## -- Get density-weighted centroids -- ##
weighted_centroids <- dat_clean %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  group_by(IUCN) %>%
  group_modify(~ {
    kde <- MASS::kde2d(
      .x$all_aln_FROH_percent,
      .x$NROH_normalized,
      n = 100
    )
    x_idx <- findInterval(.x$all_aln_FROH_percent, kde$x)
    y_idx <- findInterval(.x$NROH_normalized, kde$y)
    x_idx <- pmax(1, pmin(x_idx, length(kde$x) - 1))
    y_idx <- pmax(1, pmin(y_idx, length(kde$y) - 1))
    weights <- kde$z[cbind(x_idx, y_idx)]
    weights <- weights / sum(weights)
    data.frame(
      cx = sum(.x$all_aln_FROH_percent * weights),
      cy = sum(.x$NROH_normalized * weights)
    )
  }) %>%
  ungroup()

#### ---- Plot for all data ---- ####
p_overlay <- ggplot(dat_clean, aes(x = all_aln_FROH_percent, y = NROH_normalized,
                                    color = IUCN, fill = IUCN)) +
  annotate("segment",
           x = 0, xend = x_max_plot,
           y = 0, yend = y_max_ref,
           linetype = "dashed", color = "grey30",
           linewidth = 0.7, alpha = 0.8) +
  geom_polygon(
    data = roh_hulls,
    aes(fill = IUCN, color = IUCN),
    alpha = 0.10,
    linewidth = 0.7
  ) +
  geom_point(
    data = weighted_centroids,
    aes(x = cx, y = cy, fill = IUCN),
    size = 6,
    shape = 21,
    color = "black",
    stroke = 1.2,
    alpha = 0.75,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = iucn_colors, guide = "none") +
  scale_fill_manual(
    values = iucn_colors,
    name = "IUCN status",
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        alpha = 0.75,
        stroke = 1.2,
        color = "black"
      )
    )
  ) +
  scale_x_sqrt(
    limits = c(x_min_plot, x_max_plot),
    expand = c(0, 0),
    breaks = c(0, 5, 10, 20, 30, 50, 75, 100),
    labels = c("0", "5", "10", "20", "30", "50", "75", "100")
  ) +
  scale_y_sqrt(
    limits = c(y_min_plot, y_max_plot),
    expand = c(0, 0),
    breaks = c(0, 100, 500, 1000, 2000, 4000, 6000, 8000),
    labels = c("0", "100", "500", "1000", "2000", "4000", "6000", "8000")
  ) +
  labs(
    x = expression(bold(F[ROH] ~ "(%)" ~ (sqrt ~ scale))),
    y = expression(bold(N[ROH]))
  ) +
  My_Theme + 
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey60",
                                     linewidth = 0.4),
    legend.key.size = unit(0.55, "cm"),
    legend.text = element_text(size = 9, color = "grey10", family = "Arial"),
    legend.title = element_text(size = 10, face = "bold", color = "grey10",
                                family = "Arial"),
    plot.margin = margin(0, 0, 5, 5)
  )

p_x <- ggplot(dat_clean, aes(x = all_aln_FROH_percent, y = IUCN,
                              fill = IUCN, color = IUCN)) +
  geom_jitter(aes(color = IUCN),
              height = 0.18, size = 2.8, alpha = 0.55, shape = 21,
              stroke = 0.6) +
  geom_density_ridges(
    alpha = 0.6,
    linewidth = 0.7,
    scale = 1.3,
    quantile_lines = TRUE,
    quantiles = 2,
    color = "black",
    panel_scaling = TRUE
  ) +
  scale_fill_manual(values = iucn_colors, guide = "none") +
  scale_color_manual(values = iucn_colors, guide = "none") +
  scale_x_sqrt(
    limits = c(x_min_plot, x_max_plot),
    expand = c(0, 0),
    breaks = c(0, 5, 10, 20, 30, 50, 75, 100),
    labels = c("0", "5", "10", "20", "30", "50", "75", "100")
  ) +
  labs(x = NULL, y = "IUCN status") +
  My_Theme + 
  theme(
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "grey20", family = "Arial"),
    plot.margin = margin(5, 0, 0, 5)
  )

p_y <- ggplot(dat_clean, aes(x = NROH_normalized, y = IUCN,
                              fill = IUCN, color = IUCN)) +
  geom_jitter(aes(color = IUCN),
              height = 0.18, size = 2.8, alpha = 0.55, shape = 21,
              stroke = 0.6) +
  geom_density_ridges(
    alpha = 0.6,
    linewidth = 0.7,
    scale = 1.3,
    quantile_lines = TRUE,
    quantiles = 2,
    color = "black",
    panel_scaling = TRUE
  ) +
  scale_fill_manual(values = iucn_colors, guide = "none") +
  scale_color_manual(values = iucn_colors, guide = "none") +
  scale_x_sqrt(
    limits = c(y_min_plot, y_max_plot),
    expand = c(0, 0),
    breaks = c(0, 100, 500, 1000, 2000, 4000, 6000, 8000),
    labels = c("0", "100", "500", "1000", "2000", "4000", "6000", "8000")
  ) +
  coord_flip() +
  labs(x = expression(bold(N[ROH])), y = NULL) +
  My_Theme + 
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "grey20", family = "Arial"),
    plot.margin = margin(0, 5, 5, 0)
  )

layout <- "
AAAAAAB
CCCCCCD
"

p_combined <- p_x + plot_spacer() + p_overlay + p_y +
  plot_layout(
    design = layout,
    widths = c(5, 1),
    heights = c(1, 2.5)
  ) +
  plot_annotation(
    title = "ROH distribution across IUCN conservation categories",
    subtitle = "Convex hulls with density-weighted centroids and marginal distributions",
    theme = theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 14, face = "bold", color = "grey10",
                                hjust = 0, family = "Arial",
                                margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0,
                                   face = "italic", family = "Arial",
                                   margin = margin(b = 8)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 12, 12, 12)
    )
  )

contours_all_file_name <- paste0(results_directory, "Supplementary_Figures/", date, "_IUCN_ROH_all_density_contours.svg")
ggsave(
  contours_all_file_name, 
  p_combined,
  width = 18, 
  height = 14, 
  device = svglite)

#### ---- Run statistical tests for all data ---- ####
stats_all_output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_IUCN_ROH_all_stats.txt")

## -- Prepare scaled data -- ##
dat_ord <- dat_clean %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  mutate(
    froh_sqrt = sqrt(all_aln_FROH_percent),
    nroh_sqrt = sqrt(NROH_normalized),
    froh_scaled = scale(froh_sqrt),
    nroh_scaled = scale(nroh_sqrt)
  )

dist_mat <- dist(cbind(dat_ord$froh_scaled, dat_ord$nroh_scaled))

## -- PERMANOVA -- ##
set.seed(42)
permanova_result <- adonis2(
  dist_mat ~ IUCN,
  data = dat_ord,
  permutations = 999,
  by = "margin"
)

## -- Betadisper + permutest -- ##
betadisp_result <- betadisper(dist_mat, dat_ord$IUCN)
betadisp_test <- permutest(betadisp_result, permutations = 999)
betadisp_tukey <- TukeyHSD(betadisp_result)

## -- Pairwise PERMANOVA -- ##
iucn_levels <- levels(dat_ord$IUCN)
pairs <- combn(iucn_levels, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(pair) {
  idx <- dat_ord$IUCN %in% pair
  sub_dist <- as.dist(as.matrix(dist_mat)[idx, idx])
  sub_data <- dat_ord[idx, ]
  set.seed(42)
  res <- adonis2(sub_dist ~ IUCN, data = sub_data, permutations = 999)
  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    R2 = round(res$R2[1], 3),
    F_stat = round(res$F[1], 3),
    p_value = res$`Pr(>F)`[1]
  )
})

pairwise_df <- do.call(rbind, pairwise_results) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "bonferroni"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

## -- Convex hull areas -- ##
hull_areas <- dat_clean %>%
  filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
  filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
  group_by(IUCN) %>%
  filter(n() >= 3) %>%
  summarise(
    n = n(),
    hull_area_raw = {
      h <- chull(all_aln_FROH_percent, NROH_normalized)
      pts <- cbind(all_aln_FROH_percent[h], NROH_normalized[h])
      x <- pts[, 1]
      y <- pts[, 2]
      abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
    },
    hull_area_sqrt = {
      h <- chull(sqrt(all_aln_FROH_percent), sqrt(NROH_normalized))
      pts <- cbind(sqrt(all_aln_FROH_percent)[h], sqrt(NROH_normalized)[h])
      x <- pts[, 1]
      y <- pts[, 2]
      abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
    },
    .groups = "drop"
  ) %>%
  arrange(desc(hull_area_raw))

## -- Write all results to file -- ##
sink(stats_all_output_file)
cat("=== PERMANOVA (centroid differences) ===\n")
print(permanova_result)
cat("\n=== Betadisper (dispersion/hull size differences) ===\n")
print(betadisp_test)
cat("\n=== Pairwise Tukey HSD on dispersion ===\n")
print(betadisp_tukey)
cat("\n=== Pairwise PERMANOVA (Bonferroni corrected) ===\n")
print(pairwise_df)
cat("\n=== Convex hull areas per IUCN group ===\n")
print(hull_areas)
sink()

cat("Stats written to:", stats_all_output_file, "\n")
#### ---- END ---- ####




#### ---- VGP NROH/SROH plots for wild data, by lineage ---- ####
## -- Get unique lineages -- ##
lineages <- unique(dat_clean_wild$Extended_lineage)
lineages <- lineages[!is.na(lineages)]

## -- Refined base theme -- ##
My_Theme_Refined <- My_Theme +
  theme(
    panel.background   = element_rect(fill = "#FAFAFA", color = NA),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "grey40", linewidth = 0.5),
    axis.ticks         = element_line(color = "grey40", linewidth = 0.4),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.title         = element_text(size = 9, face = "bold", color = "grey20", family = "Arial"),
    axis.text          = element_text(size = 8, color = "grey30", family = "Arial"),
    plot.background    = element_rect(fill = "white", color = NA)
  )

## -- Build one plot per lineage -- ##
lineage_plots <- lapply(lineages, function(lin) {

  dat_lin <- dat_clean_wild %>% filter(Extended_lineage == lin)

  ## -- Per-lineage axis limits -- ##
  x_max <- max(dat_lin$all_aln_FROH_percent, na.rm = TRUE)
  y_max <- max(dat_lin$NROH_normalized,         na.rm = TRUE)

  x_buffer <- x_max * 0.05
  y_buffer <- y_max * 0.05

  x_min_plot <- 0
  x_max_plot <- x_max + x_buffer
  y_min_plot <- 0
  y_max_plot <- y_max + y_buffer

  ref_x_end   <- x_max_plot
  ref_y_end   <- x_max_plot
  ref_label_x <- ref_x_end * 0.70
  ref_label_y <- ref_x_end * 0.70

  ## -- Convex hulls -- ##
  roh_hulls <- dat_lin %>%
    group_by(IUCN) %>%
    filter(n() >= 3) %>%
    filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
    filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
    slice(chull(all_aln_FROH_percent, NROH_normalized))

  ## -- Density-weighted centroids -- ##
  weighted_centroids <- dat_lin %>%
    filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
    filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized)) %>%
    group_by(IUCN) %>%
    filter(n() >= 2) %>%
    group_modify(~ {
      .x <- .x %>%
        filter(!is.na(all_aln_FROH_percent) & !is.na(NROH_normalized)) %>%
        filter(is.finite(all_aln_FROH_percent) & is.finite(NROH_normalized))
      if (nrow(.x) < 2) return(data.frame(cx = NA_real_, cy = NA_real_))
      kde   <- MASS::kde2d(.x$all_aln_FROH_percent, .x$NROH_normalized, n = 100)
      x_idx <- findInterval(.x$all_aln_FROH_percent, kde$x)
      y_idx <- findInterval(.x$NROH_normalized,         kde$y)
      x_idx <- pmax(1, pmin(x_idx, length(kde$x) - 1))
      y_idx <- pmax(1, pmin(y_idx, length(kde$y) - 1))
      weights <- kde$z[cbind(x_idx, y_idx)]
      weights <- weights / sum(weights)
      data.frame(
        cx = sum(.x$all_aln_FROH_percent * weights),
        cy = sum(.x$NROH_normalized         * weights)
      )
    }) %>%
    filter(!is.na(cx)) %>%
    ungroup()

  ## -- Shared scale definitions (reused across p_overlay and p_x) -- ##
  x_scale <- scale_x_sqrt(
    limits = c(x_min_plot, x_max_plot),
    expand = c(0, 0),
    breaks = c(0, 5, 10, 20, 30, 50, 75, 100),
    labels = c("0", "5", "10", "20", "30", "50", "75", "100")
  )

  y_scale <- scale_y_sqrt(
    limits = c(y_min_plot, y_max_plot),
    expand = c(0, 0),
    breaks = c(0, 100, 500, 1000, 2000, 4000, 6000, 8000),
    labels = c("0", "100", "500", "1k", "2k", "4k", "6k", "8k")
  )

  y_scale_x <- scale_x_sqrt(   
    limits = c(y_min_plot, y_max_plot),
    expand = c(0, 0),
    breaks = c(0, 100, 500, 1000, 2000, 4000, 6000, 8000),
    labels = c("0", "100", "500", "1k", "2k", "4k", "6k", "8k")
  )

  ## -- Main scatter + hull overlay -- ##
  p_overlay <- ggplot(dat_lin,
                      aes(x = all_aln_FROH_percent, y = NROH_normalized,
                          color = IUCN, fill = IUCN)) +
    annotate("segment",
             x = 0, xend = ref_x_end,
             y = 0, yend = ref_y_end,
             linetype = "dashed", color = "#2C5F8A",
             linewidth = 0.75, alpha = 0.85) +
    annotate("label",
             x = ref_label_x, y = ref_label_y,
             label = "slope~'= 1'~(F[ROH]~'='~N[ROH])",
             parse = TRUE,
             size = 2.6, family = "Arial",
             color = "#2C5F8A", fill = "white",
             linewidth = 0.3, alpha = 0.9) +
    geom_polygon(
      data      = roh_hulls,
      aes(fill = IUCN, color = IUCN),
      alpha     = 0.12,
      linewidth = 0.65
    ) +
    geom_point(
      aes(color = IUCN),
      size = 1.8, shape = 16, alpha = 0.35
    ) +
    geom_point(
      data        = weighted_centroids,
      aes(x = cx, y = cy, fill = IUCN),
      size        = 3.5, shape = 21, color = "white",
      stroke      = 1.5, alpha = 0.92,
      inherit.aes = FALSE
    ) +
    geom_point(
      data        = weighted_centroids,
      aes(x = cx, y = cy, color = IUCN),
      size        = 3.5, shape = 21, fill = NA,
      stroke      = 0.8, alpha = 1,
      inherit.aes = FALSE
    ) +
    annotate("text",
             x = x_max_plot * 0.97, y = y_max_plot * 0.03,
             label = lin,
             hjust = 1, vjust = 0,
             size = 3.2, fontface = "bold.italic",
             color = "grey25", family = "Arial") +

    scale_color_manual(values = iucn_colors, guide = "none") +
    scale_fill_manual(values  = iucn_colors, guide = "none") +
    x_scale +
    y_scale +
    labs(
      x = expression(bold(F[ROH] ~ "(%)")),
      y = expression(bold(N[ROH]))
    ) +
    My_Theme_Refined +
    theme(
      legend.position = "none",
      plot.margin     = margin(0, 0, 5, 5) 
    )

  ## -- Marginal: FROH -- ##
  p_x <- ggplot(dat_lin,
                aes(x = all_aln_FROH_percent, y = IUCN,
                    fill = IUCN, color = IUCN)) +
    geom_jitter(aes(color = IUCN),
                height = 0.18, size = 2.2, alpha = 0.45, shape = 21,
                stroke = 0.5) +
    geom_density_ridges(
      alpha          = 0.55, linewidth = 0.6, scale = 1.2,
      quantile_lines = TRUE, quantiles = 2,
      color          = "grey20", panel_scaling = TRUE
    ) +
    scale_fill_manual(values  = iucn_colors, guide = "none") +
    scale_color_manual(values = iucn_colors, guide = "none") +
    x_scale +                               
    labs(x = NULL, y = expression(italic("IUCN"))) +
    My_Theme_Refined +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.line.x        = element_blank(),
      axis.text.y        = element_text(size = 7.5, color = "grey25", family = "Arial"),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(5, 0, 0, 5)  
    )

  ## -- Marginal: NROH -- ##
  p_y <- ggplot(dat_lin,
                aes(x = NROH_normalized, y = IUCN,
                    fill = IUCN, color = IUCN)) +
    geom_jitter(aes(color = IUCN),
                height = 0.18, size = 2.2, alpha = 0.45, shape = 21,
                stroke = 0.5) +
    geom_density_ridges(
      alpha          = 0.55, linewidth = 0.6, scale = 1.2,
      quantile_lines = TRUE, quantiles = 2,
      color          = "grey20", panel_scaling = TRUE
    ) +
    scale_fill_manual(values  = iucn_colors, guide = "none") +
    scale_color_manual(values = iucn_colors, guide = "none") +
    y_scale_x +                       
    coord_flip() +
    labs(x = expression(bold(N[ROH])), y = NULL) +
    My_Theme_Refined +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.line.y        = element_blank(),
      axis.text.x        = element_text(size = 7.5, color = "grey25", family = "Arial"),
      panel.grid.major.y = element_blank(),
      plot.margin        = margin(0, 5, 5, 0)
    )

  ## -- Assemble: use plot_layout with NULL spacer and fixed heights/widths -- ##
  layout <- "
  AAAAAAB
  CCCCCCD
  "

  p_block <- p_x + plot_spacer() + p_overlay + p_y +
    plot_layout(
      design  = layout,
      widths  = c(1, 0),    
      heights = c(1, 2.5)
    ) +
    plot_annotation(
      theme = theme(
        text            = element_text(family = "Arial"),
        plot.background = element_rect(fill = "white", color = "grey85",
                                       linewidth = 0.5),
        plot.margin     = margin(8, 8, 8, 8)
      )
    )

  return(p_block)
})

names(lineage_plots) <- lineages

## -- Shared legend -- ##
legend_source <- ggplot(dat_clean_wild,
                        aes(x = all_aln_FROH_percent, y = NROH_normalized,
                            fill = IUCN)) +
  geom_point(shape = 21, size = 4, alpha = 0.85, stroke = 1.5, color = "black") +
  scale_fill_manual(
    values = iucn_colors,
    name   = "IUCN status",
    guide  = guide_legend(
      override.aes = list(shape = 21, size = 4, alpha = 0.85,
                          stroke = 1.5, color = "black")
    )
  ) +
  My_Theme_Refined +
  theme(
    legend.position   = "right",
    legend.background = element_rect(fill = "white", color = "grey70",
                                     linewidth = 0.4),
    legend.margin     = margin(3, 6, 3, 4),
    legend.key.size   = unit(0.4, "cm"),
    legend.spacing.y  = unit(0.1, "cm"),
    legend.text       = element_text(size = 8,  color = "grey10", family = "Arial"),
    legend.title      = element_text(size = 9, face = "bold",
                                     color = "grey10", family = "Arial")
  )

shared_legend <- cowplot::get_legend(legend_source)

## -- Assemble final figure -- ##
n_lin  <- length(lineages)
ncols  <- ceiling(sqrt(n_lin))
nrows  <- ceiling(n_lin / ncols)

p_all_lineages <- wrap_plots(lineage_plots, ncol = ncols)

p_final <- (p_all_lineages | wrap_elements(shared_legend)) +
  plot_layout(widths = c(ncols * 10, 1)) +
  plot_annotation(
    title    = "ROH distribution across IUCN conservation categories by Extended Lineage",
    subtitle = "Convex hulls · density-weighted centroids · marginal distributions  |  dashed line: slope = 1",
    theme    = theme(
      text            = element_text(family = "Arial"),
      plot.title      = element_text(
        size   = 14, face = "bold", color = "grey10",
        hjust  = 0,  family = "Arial",
        margin = margin(b = 4)
      ),
      plot.subtitle   = element_text(
        size   = 9.5, color = "grey45", hjust = 0,
        face   = "italic", family = "Arial",
        margin = margin(b = 10)
      ),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(14, 14, 14, 14)
    )
  )

## -- Save Figure -- ##
panel_size <- 5
leg_width  <- 1.2

fig_w <- min(panel_size * ncols + leg_width, 49)
fig_h <- min(panel_size * nrows,             49)

contours_lineage_file_name <- paste0(results_directory, "Supplementary_Figures/", date, "_IUCN_ROH_wild_by_lineage.svg")
ggsave(
  contours_lineage_file_name,
  p_final,
  width     = fig_w,
  height    = fig_h,
  units     = "in",
  dpi       = 600,
  device    = svglite,
  limitsize = FALSE
)
#### ---- END ---- ####




#### ---- MSMC plot for all data ---- ####
## -- Prep the data -- ##
data_read <- map2_dfr(
  data$MSMC_file_name,
  data$species,
  ~{
    if (!file.exists(.x)) {
      message("Skipping missing file: ", .x)
      return(NULL)
    }
    spec_dat <- readr::read_table(.x, col_names = TRUE) %>%
  mutate(right_time_boundary = as.numeric(right_time_boundary))
    tibble(
      Generations_end   = spec_dat$left_time_boundary / mu,
      Generations_start = spec_dat$right_time_boundary / mu, 
      mid_gen           = (spec_dat$left_time_boundary + spec_dat$right_time_boundary) / (2 * mu),
      Ne                = (1 / spec_dat$lambda) / (2 * mu),
      Species = .y
    )
  }
) %>%
  left_join(dplyr::select(data, species, IUCN),
            by = c("Species" = "species"))
data_read$IUCN <- factor(data_read$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))


## -- Build common log-spaced time bins -- ##
time_grid <- exp(seq(
  log(min(data_read$Generations_start[data_read$Generations_start > 0])),
  log(max(data_read$Generations_end)),
  length.out = 200
))

## -- Bin all species -- ##
data_interp <- data_read %>%
  group_by(Species, IUCN) %>%
  group_modify(~{
    tibble(
      Generations = time_grid,
      Ne = approx(
        x    = .x$mid_gen,
        y    = .x$Ne,
        xout = time_grid,
        rule = 2
      )$y
    )
  }) %>%
  ungroup()

#### ---- Summarize data and save file ---- ####
data_summary <- data_interp %>%
  group_by(IUCN, Generations) %>%
  summarise(
    median_Ne = median(Ne, na.rm = TRUE),
    lq_Ne        = quantile(Ne, 0.25, na.rm = TRUE),
    uq_Ne        = quantile(Ne, 0.75, na.rm = TRUE),
    .groups   = "drop"
  )

write.csv(data_summary, paste0(results_directory, date, "_all_MSMC_grouped_stats.csv"))


## -- Plot -- ##
all_iucn <- ggplot(
  data_summary, 
  aes(
    x = Generations, 
    y = median_Ne,
    colour = IUCN, 
    fill = IUCN)) +
  geom_ribbon(
    aes(ymin = lq_Ne, ymax = uq_Ne), 
    alpha = 0.3, 
    colour = NA) +
  geom_line(
    linewidth = 1) +
  scale_color_manual(
    values = iucn_colors, 
    drop = FALSE
  ) + 
  scale_fill_manual(
    values = iucn_colors, 
    drop = FALSE
  ) + 
  scale_x_log10(labels = scales::label_scientific()) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(
    x      = "Generations ago",
    y      = "Effective population size (Ne)",
    colour = "IUCN status",
    fill   = "IUCN status"
  ) +
  My_Theme + 
  annotation_logticks(sides = "bl")

## -- Save file -- ##
MSMC_largex_file_name <- paste0(results_directory, "Supplementary_Figures/", date, "_MSMC_all_dat.svg")
ggsave(
  MSMC_largex_file_name,
  plot = all_iucn,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600, 
  device=svglite 
)
#### ---- END ---- ####




#### ---- MSMC contemporary Ne jitterplot, all data ---- ####
modern_dat <- data_read[which(data_read$Generations_end== 0), ]
modern_dat$IUCN <- factor(modern_dat$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))

p_modern_ne <- ggplot(
  modern_dat, 
  aes(
    x=IUCN, 
    y=Ne, 
    fill=IUCN
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha=0.5, 
    color="black", 
    lwd=1.1
  ) + 
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN), 
    stroke = 0.5,
    size=8, 
    alpha = 0.8
  ) + 
  scale_fill_manual(
    values = iucn_colors, 
    drop = FALSE
  ) +  
  guides(color = guide_legend(title = "Group")) + 
  My_Theme + 
  theme(
    legend.position = "none"
  ) + 
  scale_y_continuous(
    trans = "log10", 
    expand = expansion(mult = c(0.01, 0.03))
  ) +  
  scale_x_discrete(expand = expansion(add = 0.6)) + 
  annotation_logticks(sides = "b") + 
  coord_flip()

modern_Ne_file_name <- paste0(results_directory, "Supplementary_Figures/", date, "_MSMC_all_modern_Ne.svg")
ggsave(
  modern_Ne_file_name,
  plot = p_modern_ne,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600, 
  device=svglite 
)
#### ---- END ---- ####




#### ---- FROH for all ROH by major IUCN habitat, for wild species ---- ####
## Prep data
wild_gen_hab_dat <- wild_gen_hab_dat %>% filter(!is.na(habitats))

df_major_hab <- wild_gen_hab_dat %>%
  mutate(habitats = trimws(habitats), 
    Marine_flag = if_else(str_detect(habitats, "(^|; )(9|10|11|12)_"), "yes", "no")) %>%
  separate_rows(habitats, sep = "; ") %>% 
  separate(habitats, into = c("major_habitat", "sub_habitat"), sep = "_") %>%
  group_by(across(-sub_habitat)) %>% 
  summarise(sub_habitats = paste(unique(sub_habitat), collapse = "; "),
            .groups = "drop")

df_major_hab$major_habitat <- factor(df_major_hab$major_habitat, 
  levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", NA), 
  labels=c("Forest", "Savanna", "Shrubland", "Grassland", "Wetlands (inland)", "Rocky Areas", 
  "Caves (non-aquatic)", "Desert", "Marine Neritic", "Marine Oceanic", "Marine Deep Ocean Floor", 
  "Marine Intertidal", "Marine Coastal/Supratidal", "Artificial Terrestrial", "Artificial Aquatic", 
  "Introduced Vegetation", "Other"))

highlight_habitats <- c("Marine Neritic", "Marine Oceanic", "Marine Deep Ocean Floor", "Marine Intertidal")

habitat_colors <- setNames(
  ifelse(levels(df_major_hab$major_habitat) %in% highlight_habitats, "#1A7A8A", "#C8C5BC"),
  levels(df_major_hab$major_habitat)
)


## -- Plot -- ##
maj_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x=major_habitat, 
    y = all_aln_FROH_percent, 
    fill=major_habitat)) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA,
    alpha = 0.35,      
    lwd = 0.55,        
    width = 0.65     
  ) +
  geom_jitter(
    position = position_jitter(width = 0.18, height = 0.0),
    shape = 21,
    stroke = 0.3,       
    size = 4,        
    alpha = 0.65,
    aes(fill = major_habitat)
  ) +
  scale_fill_manual(values = habitat_colors) +
  scale_y_continuous(
    name = "Inbreeding coefficient (FROH) for ROH ≥ 100 kb (%)",
    trans = "sqrt",
    breaks = seq(0, 100, by = 10),
    labels = seq(0, 100, by = 10),
    expand = expansion(mult = c(0.02, 0.04)),
    limits = c(0, NA)
  ) +
  xlab(NULL) +            
  coord_flip() +
  theme_classic(base_size = 11) +
  My_Theme + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1))

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_major_habitats_wild_frohall.svg")
ggsave(
  output_file,
  plot = maj_hab_all_aln_FROH_percent,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device=svglite
)
#### ---- END ---- ####




#### ---- FROH for all ROH by major IUCN habitat, for all species ---- ####
maj_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab,
  aes(x = major_habitat, y = all_aln_FROH_percent, fill = major_habitat)
) +
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA,
    alpha = 0.35,      
    lwd = 0.55,        
    width = 0.65     
  ) +
  geom_jitter(
    position = position_jitter(width = 0.18, height = 0.0),
    shape = 21,
    stroke = 0.3,       
    size = 4,        
    alpha = 0.65,
    aes(fill = major_habitat)
  ) +
  scale_fill_manual(values = habitat_colors) +
  scale_y_continuous(
    name = "Inbreeding coefficient (FROH) for ROH ≥ 100 kb (%)",
    trans = "sqrt",
    breaks = seq(0, 100, by = 10),
    labels = seq(0, 100, by = 10),
    expand = expansion(mult = c(0.02, 0.04)),
    limits = c(0, NA)
  ) +
  xlab(NULL) +            
  coord_flip() +
  theme_classic(base_size = 11) +
  My_Theme + 
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1))


output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_major_habitats_all_frohall.svg")
ggsave(
  output_file,
  plot = maj_hab_all_aln_FROH_percent,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device=svglite
)
#### ---- END ---- ####




#### ---- FROH for all ROH, by marine flag for wild species ---- ####
wild_gen_hab_dat <- data_wild[c(
  "species", 
  "IUCN", 
  "Combined_captivity_flag_primary", 
  "habitats", 
  "longplus_froh", 
  "all_aln_FROH_percent", 
  "heterozygosity"
  )]

wild_gen_hab_dat <- wild_gen_hab_dat %>%
  mutate(
    major_habitat_count = case_when(
      is.na(habitats) ~ NA_integer_,
      TRUE ~ habitats %>%
        str_split(";\\s*") %>%
        lapply(\(x) str_extract(x, "^[^_]+")) %>%
        sapply(\(x) length(unique(x))) %>%
        as.integer()
    ),
    all_habitat_count = case_when(
      is.na(habitats) ~ NA_integer_,
      TRUE ~ lengths(str_split(habitats, ";\\s*"))
    )
  )

#### ---- Plot species which occur in marine environments versus those that don't ---- ####
df_major_hab <- wild_gen_hab_dat %>%
  mutate(habitats = trimws(habitats),
         Marine_flag = if_else(str_detect(habitats, "(^|; )(9|10|11|12)_"), "yes", "no"))

df_major_hab <- df_major_hab %>% filter(!is.na(Marine_flag))

marine_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x = Marine_flag, 
    y = all_aln_FROH_percent, 
    fill = Marine_flag)) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.5, 
    lwd = 1.1) + 
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 0.5,
    size=8, 
    alpha = 0.8
  ) + 
  scale_fill_manual(values=Marine_flag_colors
  ) +
  coord_flip() + 
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥100kb", 
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  My_Theme + 
  theme(legend.position = "none")

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_marine_compare_wild_frohall.svg")
ggsave(
  output_file,
  plot = marine_hab_all_aln_FROH_percent,
  width = 20, 
  height = 7,
  units = "in",
  dpi = 300, 
  device=svglite
)
#### ---- END ---- ####



#### ---- Prep marine flag data for all species ---- ####
data <- rename(data, c(longplus_froh=`1%+_aln_FROH_percent`))

gen_hab_dat <- data[c(
  "species", 
  "IUCN", 
  "Combined_captivity_flag_primary", 
  "habitats", 
  "longplus_froh", 
  "all_aln_FROH_percent", 
  "heterozygosity"
  )]

gen_hab_dat <- gen_hab_dat %>%
  mutate(
    major_habitat_count = case_when(
      is.na(habitats) ~ NA_integer_,
      TRUE ~ habitats %>%
        str_split(";\\s*") %>%
        lapply(\(x) str_extract(x, "^[^_]+")) %>%
        sapply(\(x) length(unique(x))) %>%
        as.integer()
    ),
    all_habitat_count = case_when(
      is.na(habitats) ~ NA_integer_,
      TRUE ~ lengths(str_split(habitats, ";\\s*"))
    )
  )

df_major_hab <- gen_hab_dat %>%
  mutate(habitats = trimws(habitats),
         Marine_flag = if_else(str_detect(habitats, "(^|; )(9|10|11|12)_"), "yes", "no"))

df_major_hab <- df_major_hab %>% filter(!is.na(Marine_flag))
#### ---- END ---- ####


#### ---- FROH for long ROH, by marine flag, for all species ---- ####
marine_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x=Marine_flag, 
    y = longplus_froh, 
    fill=Marine_flag)) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.5, 
    lwd = 1.1) + 
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 0.5,
    size=8, 
    alpha = 0.8
  ) + 
  scale_fill_manual(values=Marine_flag_colors
  ) +
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  coord_flip() + 
  My_Theme + 
  theme(legend.position = "none")

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_marine_compare_all_frohlongplus.svg")
ggsave(
  output_file,
  plot = marine_hab_all_aln_FROH_percent,
  width = 20, 
  height = 7,
  units = "in",
  dpi = 300, 
  device=svglite
)
#### ---- END ---- ####




#### ---- FROH for all ROH, by marine flag, for all species ---- ####
marine_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x = Marine_flag,  
    y = all_aln_FROH_percent, 
    fill = Marine_flag)) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.5, 
    lwd = 1.1) + 
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 0.5,
    size=8, 
    alpha = 0.8
  ) + 
  scale_fill_manual(values=Marine_flag_colors
  ) +
  coord_flip() + 
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥100kb", 
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, NA)
  ) + 
  My_Theme + 
  theme(legend.position="none")

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_marine_compare_all_frohall.svg")
ggsave(
  output_file,
  plot = marine_hab_all_aln_FROH_percent,
  width = 20, 
  height = 7,
  units = "in",
  dpi = 300, 
  device=svglite
)
#### ---- END ---- ####





#### ---- FROH for all ROH by microhabitat with wild data ---- ####
## Set order for microhabitats in plots
data_subset_long_wild$group <- factor(data_subset_long_wild$group, levels=c("Aer", "Arb", "Fos", "Ter", "Aqu"))

microhabitat_plot <- ggplot(
  data_subset_long_wild, 
  aes(
    x = group, 
    y = all_aln_FROH_percent
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.6,
    color = "black", 
    lwd = 1.1, 
    aes(fill=group)
    ) + 
  geom_jitter(
    alpha = 0.5, 
    position = position_jitter(
      width = 0.1, 
      height = 0), 
    size = 8, 
    shape = 21, 
    color = "black", 
    aes(fill=group)
  ) +
  coord_flip(clip = "off") + 
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by = 10), 
    labels = as.character(seq(0, 100, by = 10)),
    expand = expansion(mult = c(0.01, 0.01)), 
    limits = c(0, 100)
  ) + 
  scale_fill_manual(
    values=microhabitat_colors
  ) + 
  My_Theme + 
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), 
  legend.position = "none")

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_wild_microhabitat_all_aln_FROH_percent.svg")
ggsave(
  output_file,
  plot = microhabitat_plot,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
#### ---- END ---- ####



#### ---- Prepare data for microhabitat analysis ---- ####
data_subset <- data[c("species", "IUCN", "longplus_froh", "all_aln_FROH_percent", "heterozygosity", "Fos", "Aer", "Aqu", "Ter", "Arb")]
data_subset_long <- data_subset %>%
  pivot_longer(
    cols = 6:10, 
    names_to = "group", 
    values_to = "presence"
  )
data_subset_long <- data_subset_long %>% filter(presence == 1)

## Set order for microhabitats in plots
data_subset_long$group <- factor(data_subset_long$group, levels=c("Aer", "Arb", "Fos", "Ter", "Aqu"))
#### ---- END --- ####




#### ---- FROH for long ROH by microhabitat with all data ---- ####
microhabitat_plot_longplus <- ggplot(
  data_subset_long, 
  aes(
    x = group, 
    y = longplus_froh
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.6,
    color = "black", 
    lwd = 1.1, 
    aes(fill=group)
    ) + 
  geom_jitter(
    alpha = 0.5, 
    position = position_jitter(
      width = 0.1, 
      height = 0), 
    size = 8, 
    shape = 21, 
    color = "black", 
    aes(fill=group)
  ) +
  coord_flip(clip = "off") + 
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by = 10), 
    labels = as.character(seq(0, 100, by = 10)),
    expand = expansion(mult = c(0.01, 0.01)), 
    limits = c(0, 100)
  ) + 
  scale_fill_manual(
    values=microhabitat_colors
  ) + 
  My_Theme + 
  theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), 
    legend.position = "none")

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_microhabitat_froh_longplus.svg")
ggsave(
  output_file,
  plot = microhabitat_plot_longplus,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
#### ---- END ---- ####




#### ---- FROH for all ROH by microhabitat with all data ---- ####
microhabitat_plot <- ggplot(
  data_subset_long, 
  aes(
    x = group, 
    y = all_aln_FROH_percent
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.6,
    color = "black", 
    lwd = 1.1, 
    aes(fill=group)
    ) + 
  geom_jitter(
    alpha = 0.5, 
    position = position_jitter(
      width = 0.1, 
      height = 0), 
    size = 8, 
    shape = 21, 
    color = "black", 
    aes(fill=group)
  ) +
  coord_flip(clip = "off") + 
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by = 10), 
    labels = as.character(seq(0, 100, by = 10)),
    expand = expansion(mult = c(0.01, 0.01)), 
    limits = c(0, 100)
  ) + 
  scale_fill_manual(
    values=microhabitat_colors
  ) + 
  My_Theme + 
  theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), 
    legend.position = "none")

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_microhabitat_all_aln_FROH_percent.svg")
ggsave(
  output_file,
  plot = microhabitat_plot,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
#### ---- END ---- ####




#### ---- Heterozygosity by microhabitat for all data ---- ####
microhabitat_het <- ggplot(
  data_subset_long, 
  aes(
    x = group, 
    y = heterozygosity
  )) + 
  geom_boxplot(
    stat = "summary",
    fun.data = function(x) {
        q1  <- quantile(x, 0.25)
        q3  <- quantile(x, 0.75)
        iqr <- q3 - q1
        data.frame(
            ymin   = max(min(x), q1 - 1.5 * iqr),
            lower  = q1,
            middle = mean(x),
            upper  = q3,
            ymax   = min(max(x), q3 + 1.5 * iqr)
        )
    }, 
    outlier.shape = NA, 
    alpha = 0.6,
    color = "black", 
    lwd = 1.1, 
    aes(fill=group)
    ) + 
  geom_jitter(
    alpha = 0.5, 
    position = position_jitter(
      width = 0.1, 
      height = 0), 
    size = 8, 
    shape = 21, 
    color = "black", 
    aes(fill=group)
  ) +
  coord_flip(clip = "off") + 
  scale_y_continuous(
    trans = "sqrt", 
    breaks = seq(0, 30, by=2),
    expand = expansion(mult = c(0.01, 0.01)), 
    limits = c(0, 35)
  ) + 
  scale_fill_manual(
    values=microhabitat_colors
  ) + 
  My_Theme + 
  theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), 
    legend.position = "none")

output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_all_microhabitat_het.svg")
ggsave(
  output_file,
  plot = microhabitat_het, 
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
#### ---- END ---- ####