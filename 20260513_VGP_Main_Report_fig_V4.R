#### ==== Script for plotting figures for Gardiner et al. 2026 VGP ROH analysis ==== ####
## Date: 20260513 (May 13th, 2026)
## Author: Amanda Gardiner
## Version: 4
## GOAL: Create main and supplementary figures for publication
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
library(ggrepel)
library(ggtext)
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
library(htmlwidgets)


########## ========== Load in data and variables ========== ##########
args <- commandArgs()
date <- Sys.Date()
data <- read_csv(args[6])
results_directory <- args[7]
data <- rename(data, c(longplus_froh=`1%+_aln_FROH_percent`))

########## ========== Filtering criteria ========== ##########
## -- Filter out DD/NE species for plotting and factor the remaining species -- ##
data <- data %>% filter(IUCN != "DD/NE", !is.na(IUCN))
iucn_levels <- c("CR", "EN", "VU", "NT", "LC")
data$IUCN <- factor(data$IUCN, levels=iucn_levels)

## -- Filter out invertebrates -- ##
data <- data %>% filter(data$Extended_lineage != "Cyclostomes" & data$Extended_lineage != "Other Deuterostomes" )

## -- Normalize NROH by the number of chromosomes -- ##
data <- data %>% mutate(NROH_normalized = all_aln_NROH / num_auto_chromosomes)

## -- Set categories for IDRisk -- ##
data$IDRisk_category <- cut(
  data$IDRisk_longplus,
  breaks = c(-Inf, 0.05, 0.25, 0.5, Inf),
  labels = c("Low risk", "Moderate risk", "High risk", "Extreme risk"),
  right  = TRUE
)

## -- Filter out for wild only species -- ##
data_wild <- data %>% filter(Combined_captivity_flag_primary == "no" | Combined_captivity_flag_primary == "no likely")


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

#### ---- Set IDRisk color palette ---- ####
idrisk_colors <- c(
  "Low risk"      = "#FFD700",
  "Moderate risk" = "#DAA520",
  "High risk"     = "#B8600A",
  "Extreme risk"  = "#8B1A1A"
)

#### ---- Set IUCN shapes and linetypes ---- ####
iucn_shapes <- c(
  "LC" = 25, 
  "NT" = 24, "VU" = 23, 
  "EN" = 22, "CR" = 21
  )

iucn_linetypes <- c(
  "CR" = "solid",
  "EN" = "dotdash",
  "VU" = "dashed",
  "NT" = "longdash",
  "LC" = "dotted"
)

########## ========== MAIN FIGURE 1 ========== ##########
#### ---- Heterozygosity for wild species ---- ####
het_wild <- ggplot(data_wild, aes(x=IUCN, y=heterozygosity, color=IUCN, fill=IUCN)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
  ) +
  geom_jitter(
    position = position_jitter(0.15),
    shape = 21,
    color = "black",
    aes(fill = IUCN),
    stroke = 1,
    size = 8,
    alpha = 0.8
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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
    trans = "sqrt",
    breaks = seq(0, 30, by = 2),
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
    y=longplus_froh, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

#### ---- Save Figure 1 ---- ####
Fig_one_name <- paste0(results_directory, "Main_Figures/", date, "_Main_Figure_1.svg") 
ggsave(
  Fig_one_name,
  plot = Fig_one,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
########## ========== END FIGURE 1 ========== ##########




########## ========== MAIN FIGURE 2 ========== ##########
## -- Set threshold of only plotting species with FROH_all ≥5% -- ##
subset_dat_fivepercent_wild <- data_wild %>%
  filter(all_aln_FROH_percent >= 5)


## -- Calculate equations for plot -- ##
equation_labels <- subset_dat_fivepercent_wild %>%
  mutate(
    sqrt_FROH = sqrt(all_aln_FROH_percent),
    sqrt_NROH = sqrt(NROH_normalized)
  ) %>%
  group_by(IUCN) %>%
  group_map(~ {
    fit <- lm(sqrt_NROH ~ sqrt_FROH, data = .x)
    coefs <- coef(fit)
    r2    <- summary(fit)$r.squared
    tibble(
      IUCN      = .y$IUCN,
      intercept = coefs[1],
      slope     = coefs[2],
      r2        = r2
    )
  }) %>%
  bind_rows() %>%
  mutate(
    label = sprintf(
      "y = %.2fx %s %.2f  (R² = %.2f)",
      slope,
      ifelse(intercept >= 0, "+", "-"),
      abs(intercept),
      r2
    )
  )

equation_labels <- equation_labels %>%
  mutate(
    label = sprintf(
      "*y* = %.2f*x* %s %.2f &nbsp; (*R*² = %.2f)",
      slope,
      ifelse(intercept >= 0, "+", "-"),
      abs(intercept),
      r2
    )
  )

label_positions <- subset_dat_fivepercent_wild %>%
  group_by(IUCN) %>%
  summarise(
    x_pos = max(all_aln_FROH_percent, na.rm = TRUE) * 0.9,
    y_pos = max(NROH_normalized, na.rm = TRUE) * 0.95,
    label = first(IUCN)   
  )

## -- Plot -- ##
p_scatter <- ggplot(
  subset_dat_fivepercent_wild, 
  aes(
    x = all_aln_FROH_percent, 
    y = NROH_normalized,
    fill = IUCN, 
    shape = IUCN)
  ) +
  geom_point(
    size = 6,
    stroke = 1,
    alpha = 0.7, 
    color = "black"
  ) +
  geom_smooth(
      aes(color = IUCN, fill = IUCN, linetype = IUCN),
      method = "lm",
      se = TRUE,       
      linewidth = 1.5,
      alpha = 0.15     
    ) +
  scale_fill_manual(
    values = iucn_colors, 
    drop = FALSE
  ) +  
  scale_color_manual(
    values = iucn_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = iucn_shapes,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = iucn_linetypes,
    drop = FALSE
  ) +
  guides(color = guide_legend(title = "Group")) + 
  scale_x_sqrt(
    limits = c(3, NA), 
    expand = expansion(mult = c(0.01, 0.03)), 
    breaks = seq(0, 100, by=10), 
    labels = c("5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"),
    name = "Inbreeding Coefficient (%) for ROH ≥100kb of their chromosome"
  ) +
  scale_y_sqrt(
    limits = c(0, NA),
    expand = expansion(mult = c(0.01, 0.03)),
    breaks = c(0, 25, 50, 100, 200, 300, 400, 500, 600),
    labels = c("0", "25", "50", "100", "200", "300", "400", "500", "600"), 
    name = "Number of ROH (NROH)"
  ) +
  labs(
    x = expression(bold(F[ROH] ~ "(%)" ~ (sqrt ~ scale))),
    y = expression(bold(N[ROH]))
  ) +
  My_Theme

## -- Save figure -- ##
scatterplot_figname <- paste0(results_directory, "Main_Figures/", date, "_Main_Figure_2.svg") 
ggsave(
  scatterplot_figname,
  plot = p_scatter,
  width = 10, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
########## ========== END FIGURE 2========== ##########




########## ========== MAIN FIGURE 3 ========== ##########
#### ---- Set variables for MSMC ---- ####
mu <- 1.25e-8
cutoff_old <- 3000
cutoff_recent <- 1000

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


## -- Interpolate all species onto the bins -- ##
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
  geom_vline(xintercept = cutoff_old, linetype = "dashed") + 
  geom_vline(xintercept = cutoff_recent, linetype = "dashed") + 
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
#### ---- End ---- ####

#### ---- Jitterplot of contemoprary NE across IUCN status ---- ####
modern_dat <- data_read %>%
  group_by(Species, IUCN) %>%
  slice_min(Generations_end, n = 1) %>%   
  ungroup()
modern_dat$IUCN <- factor(modern_dat$IUCN, levels = c("CR", "EN", "VU", "NT", "LC"))

modern_dat$IUCN <- factor(modern_dat$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))

p_modern_ne <- ggplot(
  modern_dat, 
  aes(
    x=IUCN, 
    y=Ne, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

#### ---- Jitterplot of older NE across IUCN status ---- ####
old_dat <- data_read %>%
  group_by(Species, IUCN) %>%
  slice_min(abs(mid_gen - cutoff_old), n = 1) %>% 
  ungroup()

old_dat$IUCN <- factor(old_dat$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))

p_old_ne <- ggplot(
  old_dat, 
  aes(
    x=IUCN, 
    y=Ne, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

#### ---- Assemble Figure 3 ---- ####
Fig_three <- (
    wild_iucn /
    p_modern_ne / 
    p_old_ne
)

## Save figure
Fig_three_name <- paste0(results_directory, "Main_Figures/", date, "_Main_Figure_3.svg") 
ggsave(
  Fig_three_name,
  plot = Fig_three,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
########## ========== END FIGURE 3 ========== ##########




########## ========== MAIN FIGURE 4 ========== ##########
#### ---- Get habitat information for species ---- ####
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

## -- Set order for microhabitats in plots -- ##
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
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
  ) +
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 1,
    size=8, 
    alpha = 0.8
  ) + 
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
  ) +
  scale_fill_manual(values=Marine_flag_colors
  ) +
  scale_y_continuous(
    name = "Inbreeding Coefficient (%) for ROH ≥1% of their chromosome", 
    trans = "sqrt", 
    breaks = seq(0, 100, by=10), 
    labels = c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), 
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, 70)
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
    outlier.shape = NA,
    alpha = 0.5,
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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
    expand = expansion(mult = c(0.01, 0.03)), 
    limits = c(0, 70)
  ) + 
  scale_fill_manual(values = microhabitat_colors) + 
  My_Theme + 
  theme(plot.margin = margin(t = 5, r = 80, b = 5, l = 5, unit = "pt"))
#### ---- END ---- ####

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
    outlier.shape = NA,
    alpha = 0.5,
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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
  ) &
  theme(
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 80, b = 5, l = 5, unit = "pt")
  )


## -- Save figure -- ##
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
########## ========== END FIGURE 4 ========== ##########




########## ========== MAIN FIGURE 5 ========== ##########
## -- Set equation for IDRisk clines -- ##
cline_values <- c(0.05, 0.25, 0.5)

cline_dat_wild <- do.call(rbind, lapply(cline_values, function(k) {
  x_seq <- seq(min(data_wild$longplus_froh/100, na.rm = TRUE),
               max(data_wild$longplus_froh/100, na.rm = TRUE),
               length.out = 200)
  data.frame(x = x_seq, y = k / x_seq, k = as.factor(k))
}))
cline_dat_wild$k <- as.numeric(as.character(cline_dat_wild$k))
cline_meta <- data.frame(
  k      = c(0.05, 0.25, 0.5),       
  color  = c("#DAA520", "#B8600A", "#8B1A1A")  
)
cline_dat_wild <- do.call(rbind, lapply(cline_values, function(k) {
  x_seq <- seq(min(data_wild$longplus_froh/100, na.rm = TRUE),
               max(data_wild$longplus_froh/100, na.rm = TRUE),
               length.out = 200)
  data.frame(x = x_seq, y = k / x_seq, k = k)
}))
cline_dat_wild <- merge(cline_dat_wild, cline_meta, by = "k")
label_x <- max(data_wild$longplus_froh/100, na.rm = TRUE) * 0.85
region_labels <- data.frame(
  x     = rep(label_x, 4),
  y     = c(
    0.5 * cline_values[1] / label_x,                             
    sqrt(cline_values[1]/label_x * cline_values[2]/label_x),      
    sqrt(cline_values[2]/label_x * cline_values[3]/label_x),  
    1.5 * cline_values[3] / label_x                              
  ),
  label = c("Low risk", "Moderate risk", "High risk", "Extreme risk"),
  color = c("#FFD700", "#DAA520", "#B8600A", "#8B1A1A")
)

## -- Define outliers which should be labeled -- ##
data_wild$outlier_IDRisk_longplus <- ifelse(
  data_wild$IDRisk_longplus >= 0.5,
  data_wild$species,
  NA
)

## -- Plot -- ##
IDRisk_wild_IUCN_filename <- paste0("Figures/", date, '_IDRisk_longplus_wild_shape_IUCN.png')
png(file = IDRisk_wild_IUCN_filename, width = 2000, height = 2000, res = 300)

ggplot(
  data_wild,
  aes(
    x     = longplus_froh / 100,
    y     = heterozygosity,
    shape = IUCN,
    color = IDRisk_category
  )
) +
  geom_line(
    data        = cline_dat_wild,
    aes(x = x, y = y, group = k),
    color       = cline_dat_wild$color,
    linewidth   = 0.8,
    inherit.aes = FALSE
  ) +
  geom_point(size = 2) +
  scale_color_manual(
    values = idrisk_colors,
    name   = expression(ID[risk] ~ category),
    drop   = FALSE
  ) +
  geom_text(
    data        = region_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    color       = region_labels$color,
    hjust       = 1,
    size        = 3,
    fontface    = "bold"
  ) +
  xlab('FROH for ROH ≥1% of their chromosome') +
  ylab('Het/kb in non-ROH regions') +
  scale_shape_manual(values = iucn_shapes) +
  scale_x_sqrt() +
  scale_y_sqrt(limits = c(0, 35)) +
  geom_text_repel(
    data         = subset(data_wild, !is.na(outlier_IDRisk_longplus)),
    aes(label    = species),
    na.rm        = TRUE,
    size         = 3,
    max.overlaps = Inf,
    color        = "black"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line       = element_blank()
  ) +
  My_Theme

dev.off()
########## ========== END FIGURE 5 ========== ##########









########## ========== ========== ##########
########## SUPPLEMENTARY FIGURES 
########## ========== ========== ##########


########## ========== SFIG1: HETEROZYGOSITY FOR ALL SPECIES ========== ##########
het_alldat <- ggplot(
  data, 
  aes(x=IUCN, y=heterozygosity, color=IUCN, fill=IUCN)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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


## -- Save figure -- ##
Het_alldat_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_1.svg") 
ggsave(
  Het_alldat_figname,
  plot = het_alldat,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300,
  device=svglite 
)
########## ========== END ========== ##########


########## ========== SFIG2: FROH FOR ROH ≥1% OF THEIR CHROM, IN ALL DATA ========== ##########
p_froh_longplus_alldat <- ggplot(
  data, 
  aes(
    x=IUCN, 
    y=longplus_froh, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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


## -- Save figure -- ## 
froh_long_alldat_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_2.svg") 
ggsave(
  froh_long_alldat_figname,
  plot = p_froh_longplus_alldat,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300,
  device=svglite 
)
########## ========== END ========== ##########


########## ========== SFIG3: FROH FOR ROH ≥100kb, IN ALL DATA  ========== ##########
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## -- Save -- ##
froh_all_alldat_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_3.svg") 
ggsave(
  froh_all_alldat_figname,
  plot = p_frohall_alldat,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300,
  device=svglite 
)
########## ========== END ========== ##########


########## ========== SFIG4: NROH/FROH PLOT FOR ROH≥100KB, IN ALL DATA ========== ##########
## -- Set threshold of only plotting species with FROH_all ≥5% -- ##
subset_dat_fivepercent_all <- data %>%
  filter(all_aln_FROH_percent >= 5)


## -- Calculate equations for plot -- ##
equation_labels <- subset_dat_fivepercent_all %>%
  mutate(
    sqrt_FROH = sqrt(all_aln_FROH_percent),
    sqrt_NROH = sqrt(NROH_normalized)
  ) %>%
  group_by(IUCN) %>%
  group_map(~ {
    fit <- lm(sqrt_NROH ~ sqrt_FROH, data = .x)
    coefs <- coef(fit)
    r2    <- summary(fit)$r.squared
    tibble(
      IUCN      = .y$IUCN,
      intercept = coefs[1],
      slope     = coefs[2],
      r2        = r2
    )
  }) %>%
  bind_rows() %>%
  mutate(
    label = sprintf(
      "y = %.2fx %s %.2f  (R² = %.2f)",
      slope,
      ifelse(intercept >= 0, "+", "-"),
      abs(intercept),
      r2
    )
  )

equation_labels <- equation_labels %>%
  mutate(
    label = sprintf(
      "*y* = %.2f*x* %s %.2f &nbsp; (*R*² = %.2f)",
      slope,
      ifelse(intercept >= 0, "+", "-"),
      abs(intercept),
      r2
    )
  )

## -- Plot -- ##
p_scatter <- ggplot(
  subset_dat_fivepercent_all, 
  aes(
    x = all_aln_FROH_percent, 
    y = NROH_normalized,
    color = IUCN, 
    fill = IUCN, 
    shape = IUCN)
  ) +
  geom_point(
    aes(color = IUCN), 
    size = 6,
    color = "black",
    stroke = 1,
    alpha = 0.7
  ) +
  geom_smooth(
      aes(color = IUCN, fill = IUCN, linetype = IUCN),
      method = "lm",
      se = TRUE,       
      linewidth = 1,
      alpha = 0.15     
    ) +
  scale_fill_manual(
    values = iucn_colors, 
    drop = FALSE
  ) +  
  scale_color_manual(
    values = iucn_colors,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = iucn_linetypes,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = iucn_shapes, 
    drop = FALSE
  ) + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_sqrt(
    limits = c(5, 90), 
    expand = expansion(mult = c(0.01, 0.03)), 
    breaks = seq(0, 100, by=10), 
    labels = c("5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"),
    name = "Inbreeding Coefficient (%) for ROH ≥100kb of their chromosome"
  ) +
  scale_y_sqrt(
    limits = c(0, 1000),
    expand = expansion(mult = c(0.01, 0.03)),
    breaks = c(0, 25, 50, 100, 200, 300, 400, 600, 800, 1000),
    labels = c("0", "25", "50", "100", "200", "300", "400", "600", "800", "1000")
  ) +
  labs(
    x = expression(bold(F[ROH] ~ "(%)" ~ (sqrt ~ scale))),
    y = expression(bold(N[ROH]))
  ) +
  My_Theme

## -- Save figure -- ##
scatterplot_figname <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_4.svg") 
ggsave(
  scatterplot_figname,
  plot = p_scatter,
  width = 10, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
########## ========== END ========== ##########


########## ========== SFIG5: MSMC plot for all data ========== ##########
#### ---- Set variables for MSMC ---- ####
mu <- 1.25e-8
cutoff_old <- 3000
cutoff_recent <- 1000

#### ---- Read MSMC data ---- ####
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


## -- Build common log time bins -- ##
time_grid <- exp(seq(
  log(min(data_read$Generations_start[data_read$Generations_start > 0])),
  log(max(data_read$Generations_end)),
  length.out = 200
))


## -- Interpolate all species onto the bins -- ##
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
all_iucn <- ggplot(
  data_summary,
  aes(
    x     = Generations, 
    y     = median_Ne, 
    color = IUCN, 
    fill  = IUCN)) +
  geom_vline(xintercept = cutoff_old, linetype="dashed") + 
  geom_vline(xintercept = cutoff_recent, linetype="dashed") + 
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
#### ---- End ---- ####

#### ---- Jitterplot of contemoprary NE across IUCN status ---- ####
modern_dat <- data_read[which(data_read$Generations_end == 0), ]
modern_dat$IUCN <- factor(modern_dat$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))

p_modern_ne <- ggplot(
  modern_dat, 
  aes(
    x=IUCN, 
    y=Ne, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

#### ---- Jitterplot of older NE across IUCN status ---- ####
old_dat <- data_read[data_read$Generations_end <= 3000 & data_read$Generations_start >= 3000, ]
old_dat$IUCN <- factor(old_dat$IUCN, levels=c("CR", "EN", "VU", "NT", "LC"))

p_old_ne <- ggplot(
  old_dat, 
  aes(
    x=IUCN, 
    y=Ne, 
    fill=IUCN
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

#### ---- Assemble Figure ---- ####
Fig <- (
    all_iucn /
    p_modern_ne / 
    p_old_ne
)

## -- Save figure -- ##
msmc_supp_fig <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_5.svg") 
ggsave(
  msmc_supp_fig,
  plot = Fig,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
########## ========== END ========== ##########


########## ========== SFIG6: FROH FOR ROH≥100kb, BY MAJOR IUCN HABITAT, WILD SPECIES ========== ##########
## -- Prep data -- ##
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
  levels=c(
    "1", "2", "3", 
    "4", "5", "6", 
    "7", "8", "9", 
    "10", "11", "12", 
    "13", "14", "15", 
    "16", "17", NA), 
  labels=c(
    "Forest", "Savanna", "Shrubland", 
    "Grassland", "Wetlands (inland)", "Rocky Areas", 
    "Caves (non-aquatic)", "Desert", "Marine Neritic", 
    "Marine Oceanic", "Marine Deep Ocean Floor", "Marine Intertidal", 
    "Marine Coastal/Supratidal", "Artificial Terrestrial", "Artificial Aquatic", 
    "Introduced Vegetation", "Other"))

## -- Highlight the marine habitats on the plots -- ##
highlight_habitats <- c(
  "Marine Neritic", 
  "Marine Oceanic", 
  "Marine Deep Ocean Floor", 
  "Marine Intertidal")

habitat_colors <- setNames(
  ifelse(levels(df_major_hab$major_habitat) %in% highlight_habitats, "#1A7A8A", "#C8C5BC"),
  levels(df_major_hab$major_habitat)
)

# -- Plot -- ##
maj_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x=major_habitat, 
    y = all_aln_FROH_percent, 
    fill=major_habitat)) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1, 
    aes(fill=major_habitat)
  ) +
  geom_jitter(
    position = position_jitter(width = 0.18, height = 0.0),
    shape = 21,
    stroke = 0.3,       
    size = 4,        
    alpha = 0.65,
    aes(fill = major_habitat)
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_6.svg")
ggsave(
  output_file,
  plot = maj_hab_all_aln_FROH_percent,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device=svglite
)
########## ========== END ========== ##########


########## ========== SFIG7: FROH FOR ROH≥100kb, BY MAJOR IUCN HABITAT, ALL SPECIES ========== ##########
## -- Prep data -- ##
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

gen_hab_dat <- gen_hab_dat %>% filter(!is.na(habitats))

df_major_hab <- gen_hab_dat %>%
  mutate(habitats = trimws(habitats), 
    Marine_flag = if_else(str_detect(habitats, "(^|; )(9|10|11|12)_"), "yes", "no")) %>%
  separate_rows(habitats, sep = "; ") %>% 
  separate(habitats, into = c("major_habitat", "sub_habitat"), sep = "_") %>%
  group_by(across(-sub_habitat)) %>% 
  summarise(sub_habitats = paste(unique(sub_habitat), collapse = "; "),
            .groups = "drop")

df_major_hab$major_habitat <- factor(df_major_hab$major_habitat, 
  levels=c(
    "1", "2", "3", 
    "4", "5", "6", 
    "7", "8", "9", 
    "10", "11", "12", 
    "13", "14", "15", 
    "16", "17", NA), 
  labels=c(
    "Forest", "Savanna", "Shrubland", 
    "Grassland", "Wetlands (inland)", "Rocky Areas", 
    "Caves (non-aquatic)", "Desert", "Marine Neritic", 
    "Marine Oceanic", "Marine Deep Ocean Floor", "Marine Intertidal", 
    "Marine Coastal/Supratidal", "Artificial Terrestrial", "Artificial Aquatic", 
    "Introduced Vegetation", "Other"))

## -- Highlight the marine habitats on the plots -- ##
highlight_habitats <- c(
  "Marine Neritic", 
  "Marine Oceanic", 
  "Marine Deep Ocean Floor", 
  "Marine Intertidal")

habitat_colors <- setNames(
  ifelse(levels(df_major_hab$major_habitat) %in% highlight_habitats, "#1A7A8A", "#C8C5BC"),
  levels(df_major_hab$major_habitat)
)

## -- Plot -- ##
maj_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab,
  aes(x = major_habitat, y = all_aln_FROH_percent, fill = major_habitat)
  ) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1, 
    aes(fill=major_habitat)
  ) +
  geom_jitter(
    position = position_jitter(width = 0.18, height = 0.0),
    shape = 21,
    stroke = 0.3,       
    size = 4,        
    alpha = 0.65,
    aes(fill = major_habitat)
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_7.svg")
ggsave(
  output_file,
  plot = maj_hab_all_aln_FROH_percent,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device=svglite
)
########## ========== END ========== ##########


########## ========== SFIG8: FROH FOR ROH>100kb, BY MARINE FLAG, WILD SPECIES ========== ##########
## -- Prep data -- ##
df_major_hab <- wild_gen_hab_dat %>%
  mutate(habitats = trimws(habitats),
         Marine_flag = if_else(str_detect(habitats, "(^|; )(9|10|11|12)_"), "yes", "no"))

df_major_hab <- df_major_hab %>% filter(!is.na(Marine_flag))

## -- Plot -- ##
marine_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x = Marine_flag, 
    y = all_aln_FROH_percent, 
    fill = Marine_flag)) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1, 
    aes(fill=Marine_flag)
  ) +
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 0.5,
    size=8, 
    alpha = 0.8
  ) + 
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_8.svg")
ggsave(
  output_file,
  plot = marine_hab_all_aln_FROH_percent,
  width = 20, 
  height = 7,
  units = "in",
  dpi = 300, 
  device=svglite
)
########## ========== END ========== ##########


########## ========== SFIG9: FROH FOR ROH>1% OF THEIR CHROM, BY MARINE FLAG, ALL SPECIES ========== ##########
## -- Prep data -- ##
data <- rename(data, c(longplus_froh=longplus_froh))

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

## -- Plot -- ##
marine_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x=Marine_flag, 
    y = longplus_froh, 
    fill=Marine_flag)) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1, 
    aes(fill=Marine_flag)
  ) + 
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 0.5,
    size=8, 
    alpha = 0.8
  ) + 
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
  ) +
  scale_fill_manual(values=Marine_flag_colors) +
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_9.svg")
ggsave(
  output_file,
  plot = marine_hab_all_aln_FROH_percent,
  width = 20, 
  height = 7,
  units = "in",
  dpi = 300, 
  device=svglite
)
########## ========== END ========== ##########


########## ========== SFIG10: FROH FOR ROH>100kb OF THEIR CHROM, BY MARINE FLAG, ALL SPECIES  ========== ##########
## -- Plot -- ##
marine_hab_all_aln_FROH_percent <- ggplot(
  data = df_major_hab, 
  aes(
    x = Marine_flag,  
    y = all_aln_FROH_percent, 
    fill = Marine_flag)) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    lwd = 1.1, 
    aes(fill=Marine_flag)
  ) +
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0.0), 
    shape = 21,
    color = "black",
    aes(fill = Marine_flag), 
    stroke = 0.5,
    size=8, 
    alpha = 0.8
  ) + 
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
  ) +
  scale_fill_manual(values=Marine_flag_colors) +
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_10.svg")
ggsave(
  output_file,
  plot = marine_hab_all_aln_FROH_percent,
  width = 20, 
  height = 7,
  units = "in",
  dpi = 300, 
  device=svglite
)
########## ========== END ========== ##########


########## ========== SFIG11: FROH for ROH≥100kb BY MICROHABITAT, WILD SPECIES ========== ##########
## -- Set order for microhabitats in plots -- ##
data_subset_long_wild$group <- factor(data_subset_long_wild$group, levels=c("Aer", "Arb", "Fos", "Ter", "Aqu"))

## -- Plot -- ##
microhabitat_plot <- ggplot(
  data_subset_long_wild, 
  aes(
    x = group, 
    y = all_aln_FROH_percent
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_11.svg")
ggsave(
  output_file,
  plot = microhabitat_plot,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
########## ========== END ========== ##########


########## ========== SFIG12: FROH FOR ROH≥1% OF THEIR CHROM, ALL SPECIES ========== ##########
## -- Prep data -- ##
data_subset <- data[c("species", "IUCN", "longplus_froh", "all_aln_FROH_percent", "heterozygosity", "Fos", "Aer", "Aqu", "Ter", "Arb")]
data_subset_long <- data_subset %>%
  pivot_longer(
    cols = 6:10, 
    names_to = "group", 
    values_to = "presence"
  )
data_subset_long <- data_subset_long %>% filter(presence == 1)

## -- Set order for microhabitats in plots -- ##
data_subset_long$group <- factor(data_subset_long$group, levels=c("Aer", "Arb", "Fos", "Ter", "Aqu"))

## -- Plot -- ##
microhabitat_plot_longplus <- ggplot(
  data_subset_long, 
  aes(
    x = group, 
    y = longplus_froh
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_12.svg")
ggsave(
  output_file,
  plot = microhabitat_plot_longplus,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
########## ========== END ========== ##########


########## ========== SFIG13: FROH FOR ROH≥100kb, ALL SPECIEs ========== ##########
## -- Plot -- ##
microhabitat_plot <- ggplot(
  data_subset_long, 
  aes(
    x = group, 
    y = all_aln_FROH_percent
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_13.svg")
ggsave(
  output_file,
  plot = microhabitat_plot,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
########## ========== END ========== ##########


########## ========== SFIG14: HET OURISDE ROH, ALL SPECIES ========== ##########
## -- Plot -- ##
microhabitat_het <- ggplot(
  data_subset_long, 
  aes(
    x = group, 
    y = heterozygosity
  )) + 
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
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
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23, 
    size = 5,
    fill = "white",
    color = "black",
    stroke = 1.2
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

## -- Save plot -- ##
output_file <- paste0(results_directory, "Supplementary_Figures/", date, "_Supp_Fig_14.svg")
ggsave(
  output_file,
  plot = microhabitat_het, 
  width = 20, 
  height = 10,
  units = "in",
  dpi = 300, 
  device = svglite
)
########## ========== END ========== ##########


########## ========== SFIG15: IDRISK FOR ALL SPECIES ========== ##########
## -- Set clines -- ##
cline_dat <- do.call(rbind, lapply(cline_values, function(k) {
  x_seq <- seq(min(data$longplus_froh/100, na.rm = TRUE),
               max(data$longplus_froh/100, na.rm = TRUE),
               length.out = 200)
  data.frame(x = x_seq, y = k / x_seq, k = as.factor(k))
}))

cline_dat$k <- as.numeric(as.character(cline_dat$k))
cline_meta <- data.frame(
  k      = c(0.05, 0.25, 0.5),       
  color  = c("#DAA520", "#B8600A", "#8B1A1A")  
)

label_x <- max(data$longplus_froh/100, na.rm = TRUE) * 0.85

region_labels <- data.frame(
  x     = rep(label_x, 4),
  y     = c(
    0.5 * cline_values[1] / label_x,                             
    sqrt(cline_values[1]/label_x * cline_values[2]/label_x),      
    sqrt(cline_values[2]/label_x * cline_values[3]/label_x),
    1.5 * cline_values[3] / label_x      
  ),
  label = c("Low risk", "Moderate risk", "High risk", "Extreme risk"),
  color = c("#FFD700", "#DAA520", "#B8600A", "#8B1A1A")
)

## -- Label outliers -- ##
data$outlier_IDRisk_longplus <- ifelse(
  data$IDRisk_longplus >= 0.5,
  data$species,
  NA
)

## -- Plot -- ##
sfig15 <- ggplot(
  data,
  aes(
    x     = longplus_froh / 100,
    y     = heterozygosity,
    shape = IUCN,
    color = IDRisk_category
  )
) +
  geom_line(
    data         = cline_dat,
    aes(x = x, y = y, group = k),
    color        = cline_dat_wild$color,   
    linewidth    = 0.8,
    inherit.aes  = FALSE
  ) +
  geom_point(size = 2) +
  scale_color_manual(
    values = idrisk_colors,
    name   = expression(ID[risk] ~ category),
    drop   = FALSE
  ) +
  geom_text(
    data         = region_labels,
    aes(x = x, y = y, label = label),
    inherit.aes  = FALSE,
    color        = region_labels$color,
    hjust        = 1,
    size         = 3,
    fontface     = "bold"
  ) +
  xlab('FROH for ROH ≥1% of their chromosome') +
  ylab('Het/kb in non-ROH regions') +
  scale_shape_manual(values = iucn_shapes) +
  scale_x_sqrt() +
  scale_y_sqrt(limits = c(0, 35)) +
  geom_text_repel(
    data         = subset(data_wild, !is.na(outlier_IDRisk_longplus)),
    aes(label    = species),
    na.rm        = TRUE,
    size         = 3,
    max.overlaps = Inf,
    color        = "black"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line       = element_blank()
  ) +
  My_Theme


## -- Save plot -- ##
Supp_Fig_15_file_name <- paste0("Figures/", date, '_Supp_Fig_15.svg')
ggsave(
  Supp_Fig_15_file_name,
  plot = sfig15,
  width = 20, 
  height = 10,
  units = "in",
  dpi = 600,
  device=svglite 
)
########## ========== END ========== ##########


########## ========== SFIG16:   ========== ##########
########## ========== END ========== ##########


########## ========== SFIG17:   ========== ##########
########## ========== END ========== ##########


########## ========== SFIG18:   ========== ##########
########## ========== END ========== ##########


########## ========== SFIG19:   ========== ##########
########## ========== END ========== ##########