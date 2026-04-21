#### Script to plot NROH and SROH for HGDP data ####
## Date: 20260319 (March 19th, 2026)
## Author: Amanda Gardiner
## Version: 1
## GOAL: Create plots to show NROH and SROH for HGDP data to see if populations cluster
####

#### ---- Load in necessary packages ---- ####
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggpmisc)
library(ggrepel)
library(vegan)
library(svglite)

#### ---- Load in data and variables ---- ####
args <- commandArgs()
date <- Sys.Date()
data <- read_csv(args[6])
autosome_length <- as.numeric(args[7])
output_fig_name <- args[8]
output_stats_name <- args[9]
output_medroh_fig <- args[10]
output_maxroh_fig <- args[11]

#### ---- Set themes and colors for plotting ---- ####
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

#### ---- Compile NROH and SROH in bases for each individual ---- ####
roh_summary <- data %>%
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

ggsave(
  output_fig_name,
  p,
  width  = 20,
  height = 16,
  dpi    = 600,
  device = svglite
)

#### ---- Run statistical tests for population ROH data ---- ####
## -- Prepare scaled data -- ##
dat_ord <- roh_summary_clean %>%
  filter(!is.na(ROH_length_sum) & !is.na(ROH_count)) %>%
  filter(is.finite(ROH_length_sum) & is.finite(ROH_count)) %>%
  mutate(
    sroh_sqrt   = sqrt(ROH_length_sum),
    nroh_sqrt   = sqrt(ROH_count),
    sroh_scaled = as.numeric(scale(sroh_sqrt)),
    nroh_scaled = as.numeric(scale(nroh_sqrt))
  ) %>%
  filter(is.finite(sroh_scaled) & is.finite(nroh_scaled)) %>%
  mutate(Population_elastic_ID = droplevels(Population_elastic_ID))

## Drop populations with fewer than 2 individuals (adonis2 requires >= 2 per group)
pop_counts <- dat_ord %>% count(Population_elastic_ID)
valid_pops  <- pop_counts %>% filter(n >= 2) %>% pull(Population_elastic_ID)
dat_ord     <- dat_ord %>%
  filter(Population_elastic_ID %in% valid_pops) %>%
  mutate(Population_elastic_ID = droplevels(Population_elastic_ID))

cat("Individuals retained:", nrow(dat_ord), "\n")
cat("Populations retained:", nlevels(dat_ord$Population_elastic_ID), "\n")

dist_mat <- dist(cbind(dat_ord$sroh_scaled, dat_ord$nroh_scaled))

## -- PERMANOVA -- ##
set.seed(42)
permanova_result <- adonis2(
  dist_mat ~ Population_elastic_ID,
  data         = dat_ord,
  permutations = 999,
  by           = "margin"
)

## -- Betadisper + permutest -- ##
betadisp_result <- betadisper(dist_mat, dat_ord$Population_elastic_ID)
betadisp_test   <- permutest(betadisp_result, permutations = 999)
betadisp_tukey  <- TukeyHSD(betadisp_result)

## -- Pairwise PERMANOVA -- ##
pop_levels <- levels(dat_ord$Population_elastic_ID)
pairs      <- combn(pop_levels, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(pair) {
  idx      <- dat_ord$Population_elastic_ID %in% pair
  sub_dist <- as.dist(as.matrix(dist_mat)[idx, idx])
  sub_data <- dat_ord[idx, ] %>%
    mutate(Population_elastic_ID = droplevels(Population_elastic_ID))

  set.seed(42)
  res <- adonis2(sub_dist ~ Population_elastic_ID, data = sub_data, permutations = 999)

  data.frame(
    group1  = pair[1],
    group2  = pair[2],
    R2      = round(res$R2[1], 3),
    F_stat  = round(res$F[1], 3),
    p_value = res$`Pr(>F)`[1]
  )
})

pairwise_df <- do.call(rbind, pairwise_results) %>%
  mutate(
    p_adj        = p.adjust(p_value, method = "bonferroni"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

## -- Convex hull areas -- ##
hull_areas <- roh_summary_clean %>%
  filter(!is.na(ROH_length_sum) & !is.na(ROH_count)) %>%
  filter(is.finite(ROH_length_sum) & is.finite(ROH_count)) %>%
  group_by(Population_elastic_ID) %>%
  filter(n() >= 3) %>%
  summarise(
    n              = n(),
    hull_area_raw  = {
      h   <- chull(ROH_length_sum, ROH_count)
      pts <- cbind(ROH_length_sum[h], ROH_count[h])
      x   <- pts[, 1]; y <- pts[, 2]
      abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
    },
    hull_area_sqrt = {
      h   <- chull(sqrt(ROH_length_sum), sqrt(ROH_count))
      pts <- cbind(sqrt(ROH_length_sum)[h], sqrt(ROH_count)[h])
      x   <- pts[, 1]; y <- pts[, 2]
      abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
    },
    .groups = "drop"
  ) %>%
  arrange(desc(hull_area_raw))

## -- Write all results to file -- ##
sink(output_stats_name)
cat("=== PERMANOVA (centroid differences) ===\n")
print(permanova_result)
cat("\n=== Betadisper (dispersion/hull size differences) ===\n")
print(betadisp_test)
cat("\n=== Pairwise Tukey HSD on dispersion ===\n")
print(betadisp_tukey)
cat("\n=== Pairwise PERMANOVA (Bonferroni corrected) ===\n")
print(pairwise_df)
cat("\n=== Convex hull areas per population ===\n")
print(hull_areas)
sink()
cat("Stats written to:", output_stats_name, "\n")


#### ---- Plot FROH versus median ROH size in bases ---- ####
roh_summary_clean$FROH <- roh_summary_clean$ROH_length_sum / autosome_length

## Weight = number of individuals in population (n)
pop_centroids <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  summarise(
    n               = n(),
    centroid_froh   = weighted.mean(FROH, w = rep(1, n())),  
    centroid_medroh = weighted.mean(ROH_median, w = rep(1, n())),
    sd_froh         = sd(FROH),
    sd_medroh       = sd(ROH_median),
    se_froh         = sd(FROH) / sqrt(n()),
    se_medroh       = sd(ROH_median) / sqrt(n()),
    .groups         = "drop"
  )

roh_hulls_med <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  filter(n() >= 3) %>%
  slice(chull(ROH_median, FROH))

weighted_centroids_med <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  summarise(
    FROH = weighted.mean(FROH, w = ROH_count, na.rm = TRUE),
    ROH_median = weighted.mean(ROH_median, w = ROH_count, na.rm = TRUE),
    .groups = "drop"
  )

## Order populations by centroid position (SROH then NROH)
pop_order <- pop_centroids %>%
  arrange(centroid_froh, centroid_medroh) %>%
  pull(Population_elastic_ID)

roh_summary_clean <- roh_summary_clean %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

roh_hulls_med <- roh_hulls_med %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

pop_centroids <- pop_centroids %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))


p <- ggplot(roh_summary_clean, aes(x = FROH, y = ROH_median,
                                    color = Population_elastic_ID,
                                    fill  = Population_elastic_ID)) +
  geom_point(size = 1.5, alpha = 0.4) +
  geom_polygon(data = roh_hulls_med, alpha = 0.2, linewidth = 0.3) +
  geom_errorbar(
    data = pop_centroids,
    aes(x = centroid_froh,
        ymin = centroid_medroh - sd_medroh,
        ymax = centroid_medroh + sd_medroh),
    width = 0, linewidth = 0.5, inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbarh(
    data = pop_centroids,
    aes(y    = centroid_medroh,
        xmin = centroid_froh - se_froh,
        xmax = centroid_froh + se_froh),
    height = 0, linewidth = 0.5, inherit.aes = FALSE,
    color = "black"
  ) +
  geom_point(
    data = pop_centroids,
    aes(x = centroid_froh, y = centroid_medroh, fill = Population_elastic_ID),
    shape = 21, size = 3, color = "black", stroke = 0.6, inherit.aes = FALSE
  ) +
  scale_color_manual(values = pop_colors) +
  scale_fill_manual(values  = pop_colors) +
  facet_wrap(~ Population_elastic_ID) +
  labs(
    x = "FROH (Proportion of autosome in ROH)",
    y = "Median ROH length"
  ) +
  theme_classic() +
  theme(
    strip.text    = element_text(size = 6),
    axis.text     = element_text(size = 5),
    axis.title    = element_text(size = 8),
    legend.position = "none"
  )

ggsave(
  output_medroh_fig,
  p,
  width  = 20,
  height = 16,
  dpi    = 600,
  device = svglite
)


#### ---- Plot FROH versus maximum ROH size in bases ---- ####
roh_hulls_max <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  filter(n() >= 3) %>%
  slice(chull(ROH_max_value, FROH))

weighted_centroids_max <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  summarise(
    FROH = weighted.mean(FROH, w = ROH_count, na.rm = TRUE),
    ROH_max_value = weighted.mean(ROH_max_value, w = ROH_count, na.rm = TRUE),
    .groups = "drop"
  )

pop_centroids <- roh_summary_clean %>%
  group_by(Population_elastic_ID) %>%
  summarise(
    n               = n(),
    centroid_froh   = weighted.mean(FROH, w = rep(1, n())),  
    centroid_maxroh = weighted.mean(ROH_max_value, w = rep(1, n())),
    sd_froh         = sd(FROH),
    sd_maxroh       = sd(ROH_max_value),
    se_froh         = sd(FROH) / sqrt(n()),
    se_maxroh       = sd(ROH_max_value) / sqrt(n()),
    .groups         = "drop"
  )

pop_order <- pop_centroids %>%
  arrange(centroid_froh, centroid_maxroh) %>%
  pull(Population_elastic_ID)

roh_summary_clean <- roh_summary_clean %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

roh_hulls_max <- roh_hulls_max %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

pop_centroids <- pop_centroids %>%
  mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

p <- ggplot(roh_summary_clean, aes(x = FROH, y = ROH_max_value,
                                    color = Population_elastic_ID,
                                    fill  = Population_elastic_ID)) +
  geom_point(size = 1.5, alpha = 0.4) +
  geom_polygon(data = roh_hulls_max, alpha = 0.2, linewidth = 0.3) +
  geom_errorbar(
    data = pop_centroids,
    aes(x = centroid_froh,
        ymin = centroid_maxroh - sd_maxroh,
        ymax = centroid_maxroh + sd_maxroh),
    width = 0, linewidth = 0.5, inherit.aes = FALSE,
    color = "black"
  ) +
  geom_errorbarh(
    data = pop_centroids,
    aes(y    = centroid_maxroh,
        xmin = centroid_froh - se_froh,
        xmax = centroid_froh + se_froh),
    height = 0, linewidth = 0.5, inherit.aes = FALSE,
    color = "black"
  ) +
  geom_point(
    data = pop_centroids,
    aes(x = centroid_froh, y = centroid_maxroh, fill = Population_elastic_ID),
    shape = 21, size = 3, color = "black", stroke = 0.6, inherit.aes = FALSE
  ) +
  scale_color_manual(values = pop_colors) +
  scale_fill_manual(values  = pop_colors) +
  facet_wrap(~ Population_elastic_ID) +
  labs(
    x = "FROH (Proportion of autosome in ROH)",
    y = "Maximum ROH length"
  ) +
  theme_classic() +
  theme(
    strip.text    = element_text(size = 6),
    axis.text     = element_text(size = 5),
    axis.title    = element_text(size = 8),
    legend.position = "none"
  )

ggsave(
  output_maxroh_fig,
  p,
  width  = 20,
  height = 16,
  dpi    = 600,
  device = svglite
)







# #### ---- Compile NROH and SROH in bases for each individual ---- ####
# roh_summary <- data %>%
#   group_by(IID, Population_elastic_ID) %>%
#   summarise(
#     ROH_count = n(),
#     ROH_length_sum = sum(length), 
#     .groups = "drop" 
#   )

# roh_summary_clean <- roh_summary %>%
#   filter(!is.na(Population_elastic_ID))

# ## Compute hulls
# roh_hulls <- roh_summary_clean %>%
#   group_by(Population_elastic_ID) %>%
#   filter(n() >= 3) %>%
#   slice(chull(ROH_length_sum, ROH_count))

# ## Compute centroids and order populations by centroid position
# ## First order by SROH then NROH
# pop_order <- roh_summary_clean %>%
#   group_by(Population_elastic_ID) %>%
#   summarise(
#     centroid_x = mean(ROH_length_sum),
#     centroid_y = mean(ROH_count)
#   ) %>%
#   arrange(centroid_x, centroid_y) %>%
#   pull(Population_elastic_ID)

# roh_summary_clean <- roh_summary_clean %>%
#   mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

# roh_hulls <- roh_hulls %>%
#   mutate(Population_elastic_ID = factor(Population_elastic_ID, levels = pop_order))

# ## Draw a diagonal line on the plots
# x_min <- min(roh_summary_clean$ROH_length_sum)
# x_max <- max(roh_summary_clean$ROH_length_sum)
# y_min <- min(roh_summary_clean$ROH_count)
# y_max <- max(roh_summary_clean$ROH_count)

# #### ---- PLOT ---- ####
# p <- ggplot(roh_summary_clean, aes(x = ROH_length_sum, y = ROH_count, color = Population_elastic_ID, fill = Population_elastic_ID)) +
#   annotate("segment",
#            x = x_min, xend = x_max,
#            y = y_min, yend = y_max,
#            linetype = "dashed", color = "grey40", linewidth = 0.4) +
#   geom_point(size = 1.5, alpha = 0.8) +
#   geom_polygon(data = roh_hulls, alpha = 0.2, linewidth = 0.3) +
#   scale_color_manual(values = pop_colors) +
#   scale_fill_manual(values = pop_colors) +
#   facet_wrap(~ Population_elastic_ID) +
#   labs(
#     x = "SROH (sum of ROH length)",
#     y = "NROH (number of ROH)"
#   ) +
#   theme_classic() +
#   theme(
#     strip.text = element_text(size = 6),
#     axis.text = element_text(size = 5),
#     axis.title = element_text(size = 8),
#     legend.position = "none"
#   )

# ggsave(
#   output_fig_name, 
#   p, 
#   width = 20, 
#   height = 16, 
#   dpi = 600, 
#   device = svglite
#   )

# #### ---- Do stats on the plots ---- ####
# ## -- Prepare scaled data -- ##
# dat_ord <- roh_summary_clean %>%
#   filter(!is.na(ROH_length_sum) & !is.na(ROH_count)) %>%
#   filter(is.finite(ROH_length_sum) & is.finite(ROH_count)) %>%
#   mutate(
#     sroh_sqrt   = sqrt(ROH_length_sum),
#     nroh_sqrt   = sqrt(ROH_count),
#     sroh_scaled = scale(sroh_sqrt),
#     nroh_scaled = scale(nroh_sqrt)
#   )  %>%
#   filter(is.finite(sroh_scaled) & is.finite(nroh_scaled))

# dist_mat <- dist(cbind(dat_ord$sroh_scaled, dat_ord$nroh_scaled))

# ## -- PERMANOVA -- ##
# set.seed(42)
# permanova_result <- adonis2(
#   dist_mat ~ Population_elastic_ID,
#   data         = dat_ord,
#   permutations = 999,
#   by           = "margin"
# )

# ## -- Betadisper + permutest -- ##
# betadisp_result <- betadisper(dist_mat, dat_ord$Population_elastic_ID)
# betadisp_test   <- permutest(betadisp_result, permutations = 999)
# betadisp_tukey  <- TukeyHSD(betadisp_result)

# ## -- Pairwise PERMANOVA -- ##
# pop_levels <- levels(dat_ord$Population_elastic_ID)
# pairs      <- combn(pop_levels, 2, simplify = FALSE)

# pairwise_results <- lapply(pairs, function(pair) {
#   idx      <- dat_ord$Population_elastic_ID %in% pair
#   sub_dist <- as.dist(as.matrix(dist_mat)[idx, idx])
#   sub_data <- dat_ord[idx, ]

#   set.seed(42)
#   res <- adonis2(sub_dist ~ Population_elastic_ID, data = sub_data, permutations = 999)

#   data.frame(
#     group1  = pair[1],
#     group2  = pair[2],
#     R2      = round(res$R2[1], 3),
#     F_stat  = round(res$F[1], 3),
#     p_value = res$`Pr(>F)`[1]
#   )
# })

# pairwise_df <- do.call(rbind, pairwise_results) %>%
#   mutate(
#     p_adj        = p.adjust(p_value, method = "bonferroni"),
#     significance = case_when(
#       p_adj < 0.001 ~ "***",
#       p_adj < 0.01  ~ "**",
#       p_adj < 0.05  ~ "*",
#       TRUE          ~ "ns"
#     )
#   )

# ## -- Convex hull areas -- ##
# hull_areas <- roh_summary_clean %>%
#   filter(!is.na(ROH_length_sum) & !is.na(ROH_count)) %>%
#   filter(is.finite(ROH_length_sum) & is.finite(ROH_count)) %>%
#   group_by(Population_elastic_ID) %>%
#   filter(n() >= 3) %>%
#   summarise(
#     n             = n(),
#     hull_area_raw = {
#       h   <- chull(ROH_length_sum, ROH_count)
#       pts <- cbind(ROH_length_sum[h], ROH_count[h])
#       x   <- pts[, 1]; y <- pts[, 2]
#       abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
#     },
#     hull_area_sqrt = {
#       h   <- chull(sqrt(ROH_length_sum), sqrt(ROH_count))
#       pts <- cbind(sqrt(ROH_length_sum)[h], sqrt(ROH_count)[h])
#       x   <- pts[, 1]; y <- pts[, 2]
#       abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
#     },
#     .groups = "drop"
#   ) %>%
#   arrange(desc(hull_area_raw))

# ## -- Write all results to file -- ##
# sink(output_stats_name)
# cat("=== PERMANOVA (centroid differences) ===\n")
# print(permanova_result)
# cat("\n=== Betadisper (dispersion/hull size differences) ===\n")
# print(betadisp_test)
# cat("\n=== Pairwise Tukey HSD on dispersion ===\n")
# print(betadisp_tukey)
# cat("\n=== Pairwise PERMANOVA (Bonferroni corrected) ===\n")
# print(pairwise_df)
# cat("\n=== Convex hull areas per population ===\n")
# print(hull_areas)
# sink()
# cat("Stats written to:", output_stats_name, "\n")
