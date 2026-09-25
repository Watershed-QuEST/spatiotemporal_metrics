#### Read me ####
# This script is used to calculate synchrony between sites within a watershed with a toy dataset for QuEST's planned commentary manuscript on the use of spatiotemporal metrics for assessing the spatiotemporal variance of stream chemistry across stream networks.
# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Eva Tipps, 2026-09-04 (with help building functions from Claude version 1.24012.9 (03c61d) 2026-09-04T09:58:04.000Z... reviewed and edited by E. Tipps)
# Last update (Person, Date): Eva Tipps, 2026-09-16

# Requires: 02_build_synthetic_data.R must be run first (produces data/nm_clean.csv, data/nm_field_setups.rds, and data/nm_synthetic_extended.csv), and spatiotemporal_helpers.R must be in the same folder as this script.



#### Libraries ####
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggtern)
library(ggpubr)
library(ggpmisc)
library(plotly)
library(lubridate)
library(tidyr)
library(tidyverse)
library(dataRetrieval)
library(purrr)
library(broom)
library(stats)

source("spatiotemporal_helpers.R")

#### Imports ####
data_out_dir <- "data"
plot_dir     <- "plots"

nm_toy           <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
nm_toy_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)
nm_field_setups     <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))

br_toy           <- read_csv(file.path(data_out_dir, "br_clean.csv"), show_col_types = FALSE)
br_toy_extended <- read_csv(file.path(data_out_dir, "br_synthetic_extended.csv"), show_col_types = FALSE)
br_field_setups     <- readRDS(file.path(data_out_dir, "br_field_setups.rds"))

# Monte Carlo iterations
n_iter    <- 200 #test
camp_grid <- c(2, 3, 4, 5, 7, 10, 15, 20, 30, 50)      # campaign counts to test
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

n_synthetic_months <- 60      # must match 02_build_synthetic_data.R
time_unit_days     <- 30.44   # must match 02_build_synthetic_data.R
total_days <- n_synthetic_months * time_unit_days  # span campaign times are drawn from

actual_n_campaigns_nm <- n_distinct(nm_toy$CampaignID)
actual_n_campaigns_br <- n_distinct(br_toy$CampaignID)

set.seed(42)

#### Temporal synchrony function ####
synchrony <- function(x, y) {
  n <- length(x)
  
  if (length(x) != length(y)) stop("x and y must have the same length")
  if (n < 2) stop("Need at least 2 observations")
  
  sum((x - mean(x)) * (y - mean(y))) / ((n - 1) * sd(x) * sd(y))
}

## Get all unique site pairs
#NM#
nm_sites <- nm_toy %>%
  pull(Site) %>%
  unique()

nm_site_pairs <- combn(nm_sites, 2, simplify = FALSE)

#BR#
br_sites <- br_toy %>%
  pull(Site) %>%
  unique()

br_site_pairs <- combn(br_sites, 2, simplify = FALSE)

#### Synchrony of one variable between two sites####
## NM ## 
# NPOC #
# Loop over each pair and calculate synchrony for one variable
results <- lapply(nm_site_pairs, function(pair) {
  
  paired_data <- nm_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, NPOC..mg.C.L.) %>%
    filter(!is.na(NPOC..mg.C.L.)) %>% 
    group_by(Date, Site) %>%
    slice_head(n = 1) %>% 
    summarise(NPOC..mg.C.L. = mean(NPOC..mg.C.L., na.rm = TRUE), .groups = "drop") %>%  # average duplicates
    pivot_wider(names_from = Site, values_from = NPOC..mg.C.L.) %>%
    drop_na()
  
  
  # Skip pairs with too few shared dates
  if (nrow(paired_data) < 2) {
    return(data.frame(site_x = pair[1], site_y = pair[2], synchrony = NA))
  }
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]])
  )
})

results_nm <- bind_rows(results)
results_nm <- na.omit(results_nm)

# TDN #
# Loop over each pair and calculate synchrony for one variable
nresults <- lapply(nm_site_pairs, function(pair) {
  
  paired_data <- nm_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, TDN..mg.N.L.) %>%
    filter(!is.na(TDN..mg.N.L.)) %>% 
    group_by(Date, Site) %>%
    slice_head(n = 1) %>% 
    summarise(TDN..mg.N.L. = mean(TDN..mg.N.L., na.rm = TRUE), .groups = "drop") %>%  # average duplicates
    pivot_wider(names_from = Site, values_from = TDN..mg.N.L.) %>%
    drop_na()
  
  
  # Skip pairs with too few shared dates
  if (nrow(paired_data) < 2) {
    return(data.frame(site_x = pair[1], site_y = pair[2], synchrony = NA))
  }
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]])
  )
})

nresults_nm <- bind_rows(nresults)
nresults_nm <- na.omit(nresults_nm)

## BR ##
# NPOC #
# Loop over each pair and calculate synchrony for one variable
results <- lapply(br_site_pairs, function(pair) {
  
  paired_data <- br_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, NPOC..mg.C.L.) %>%
    filter(!is.na(NPOC..mg.C.L.)) %>% 
    group_by(Date, Site) %>%
    slice_head(n = 1) %>% 
    summarise(NPOC..mg.C.L. = mean(NPOC..mg.C.L., na.rm = TRUE), .groups = "drop") %>%  # average duplicates
    pivot_wider(names_from = Site, values_from = NPOC..mg.C.L.) %>%
    drop_na()
  
  
  # Skip pairs with too few shared dates
  if (nrow(paired_data) < 2) {
    return(data.frame(site_x = pair[1], site_y = pair[2], synchrony = NA))
  }
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]])
  )
})

results_br <- bind_rows(results)
results_br <- na.omit(results_br)

# TDN #
# Loop over each pair and calculate synchrony for one variable
nresults <- lapply(br_site_pairs, function(pair) {
  
  paired_data <- br_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, TDN..mg.N.L.) %>%
    filter(!is.na(TDN..mg.N.L.)) %>% 
    group_by(Date, Site) %>%
    slice_head(n = 1) %>% 
    summarise(TDN..mg.N.L. = mean(TDN..mg.N.L., na.rm = TRUE), .groups = "drop") %>%  # average duplicates
    pivot_wider(names_from = Site, values_from = TDN..mg.N.L.) %>%
    drop_na()
  
  
  # Skip pairs with too few shared dates
  if (nrow(paired_data) < 2) {
    return(data.frame(site_x = pair[1], site_y = pair[2], synchrony = NA))
  }
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]])
  )
})

nresults_br <- bind_rows(nresults)
nresults_br <- na.omit(nresults_br)


#### Synchrony between two variables (NPOC and TDN) at one site ####
# NM #
# Loop over each site and calculate synchrony between NPOC and TDN
results <- lapply(nm_sites, function(site) {
  
  site_data <- nm_toy %>%
    filter(Site == site) %>%
    select(Date, NPOC..mg.C.L., TDN..mg.N.L.) %>%
    filter(!is.na(NPOC..mg.C.L.) & !is.na(TDN..mg.N.L.)) %>%
    group_by(Date) %>%
    summarise(
      NPOC = mean(NPOC..mg.C.L., na.rm = TRUE),
      TDN  = mean(TDN..mg.N.L.,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    drop_na()
  
  # Skip sites with too few shared dates
  if (nrow(site_data) < 2) {
    return(data.frame(Site = site, synchrony = NA))
  }
  
  data.frame(
    Site      = site,
    synchrony = synchrony(site_data$NPOC, site_data$TDN)
  )
})

cnresults_nm <- bind_rows(results) %>% na.omit()

# BR #
# Loop over each site and calculate synchrony between NPOC and TDN
results <- lapply(br_sites, function(site) {
  
  site_data <- br_toy %>%
    filter(Site == site) %>%
    select(Date, NPOC..mg.C.L., TDN..mg.N.L.) %>%
    filter(!is.na(NPOC..mg.C.L.) & !is.na(TDN..mg.N.L.)) %>%
    group_by(Date) %>%
    summarise(
      NPOC = mean(NPOC..mg.C.L., na.rm = TRUE),
      TDN  = mean(TDN..mg.N.L.,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    drop_na()
  
  # Skip sites with too few shared dates
  if (nrow(site_data) < 2) {
    return(data.frame(Site = site, synchrony = NA))
  }
  
  data.frame(
    Site      = site,
    synchrony = synchrony(site_data$NPOC, site_data$TDN)
  )
})

cnresults_br <- bind_rows(results) %>% na.omit()


#### Flag results with not enough data ####
# Mark low-confidence pairs (< 3 shared dates) as NA
# NM DOC #
results_nm_flagged <- lapply(nm_site_pairs, function(pair) {
  paired_data <- nm_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, NPOC..mg.C.L.) %>%
    filter(!is.na(NPOC..mg.C.L.)) %>%
    group_by(Date, Site) %>%
    slice_head(n = 1) %>%
    pivot_wider(names_from = Site, values_from = NPOC..mg.C.L.) %>%
    drop_na()
  
  n_shared <- nrow(paired_data)
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = if (n_shared < 3) NA else synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]]),
    n_shared  = n_shared
  )
}) %>% bind_rows()

# BR DOC #
results_br_flagged <- lapply(br_site_pairs, function(pair) {
  paired_data <- br_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, NPOC..mg.C.L.) %>%
    filter(!is.na(NPOC..mg.C.L.)) %>%
    group_by(Date, Site) %>%
    slice_head(n = 1) %>%
    pivot_wider(names_from = Site, values_from = NPOC..mg.C.L.) %>%
    drop_na()
  
  n_shared <- nrow(paired_data)
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = if (n_shared < 3) NA else synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]]),
    n_shared  = n_shared
  )
}) %>% bind_rows()

# NM TDN #
nresults_nm_flagged <- lapply(nm_site_pairs, function(pair) {
  paired_data <- nm_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, TDN..mg.N.L.) %>%
    filter(!is.na(TDN..mg.N.L.)) %>%
    group_by(Date, Site) %>%
    slice_head(n = 1) %>%
    pivot_wider(names_from = Site, values_from = TDN..mg.N.L.) %>%
    drop_na()
  
  n_shared <- nrow(paired_data)
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = if (n_shared < 3) NA else synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]]),
    n_shared  = n_shared
  )
}) %>% bind_rows()

# BR TDN #
nresults_br_flagged <- lapply(br_site_pairs, function(pair) {
  paired_data <- br_toy %>%
    filter(Site %in% pair) %>%
    select(Date, Site, TDN..mg.N.L.) %>%
    filter(!is.na(TDN..mg.N.L.)) %>%
    group_by(Date, Site) %>%
    slice_head(n = 1) %>%
    pivot_wider(names_from = Site, values_from = TDN..mg.N.L.) %>%
    drop_na()
  
  n_shared <- nrow(paired_data)
  
  data.frame(
    site_x    = pair[1],
    site_y    = pair[2],
    synchrony = if (n_shared < 3) NA else synchrony(paired_data[[pair[1]]], paired_data[[pair[2]]]),
    n_shared  = n_shared
  )
}) %>% bind_rows()

#### Heat map plots ####
## NM Heat map of DOC synchrony ##
# Site → area lookup
nmarea_lookup <- nm_toy %>%
  select(Site, Area.m2) %>%
  distinct() %>%
  group_by(Site) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Area.m2)  # sort by increasing area

# Sites ordered by watershed area
nmsites_by_area <- nmarea_lookup$Site

# Build symmetric df using area-ordered factor levels
nmarea_doc <- bind_rows(
  results_nm_flagged,
  results_nm_flagged %>% rename(site_x = site_y, site_y = site_x),
  data.frame(site_x = nmsites_by_area, site_y = nmsites_by_area, synchrony = 1)
) %>%
  mutate(
    site_x = factor(site_x, levels = nmsites_by_area),
    site_y = factor(site_y, levels = nmsites_by_area)
  )

p1 <- ggplot(nmarea_doc, aes(x = site_x, y = site_y, fill = synchrony)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(synchrony, 2)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Synchrony"
  ) +
  labs(
    title = "Upper Santa Fe Pairwise NPOC Synchrony",
    x = "Site (by increasing area)", y = "Site (by increasing area)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p1)

## BR Heat map of DOC synchrony ##
# Site → area lookup
brarea_lookup <- br_toy %>%
  select(Site, Area.m2) %>%
  distinct() %>%
  group_by(Site) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Area.m2)  # sort by increasing area

# Sites ordered by watershed area
brsites_by_area <- brarea_lookup$Site

# Build symmetric df using area-ordered factor levels
brarea_doc <- bind_rows(
  results_br_flagged,
  results_br_flagged %>% rename(site_x = site_y, site_y = site_x),
  data.frame(site_x = brsites_by_area, site_y = brsites_by_area, synchrony = 1)
) %>%
  mutate(
    site_x = factor(site_x, levels = brsites_by_area),
    site_y = factor(site_y, levels = brsites_by_area)
  ) %>%
  distinct(site_x, site_y, .keep_all=TRUE)

p2 <- ggplot(brarea_doc, aes(x = site_x, y = site_y, fill = synchrony)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(synchrony, 2)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Synchrony"
  ) +
  labs(
    title = "Brush Creek Pairwise DOC Synchrony",
    x = "Site (by increasing area)", y = "Site (by increasing area)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p2)

## NM Heat map of TDN synchrony ##

# Build symmetric df using area-ordered factor levels
nmarea_tdn <- bind_rows(
  nresults_nm_flagged,
  nresults_nm_flagged %>% rename(site_x = site_y, site_y = site_x),
  data.frame(site_x = nmsites_by_area, site_y = nmsites_by_area, synchrony = 1)
) %>%
  mutate(
    site_x = factor(site_x, levels = nmsites_by_area),
    site_y = factor(site_y, levels = nmsites_by_area)
  )

p3 <- ggplot(nmarea_tdn, aes(x = site_x, y = site_y, fill = synchrony)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(synchrony, 2)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Synchrony"
  ) +
  labs(
    title = "Upper Santa Fe Pairwise TDN Synchrony",
    x = "Site (by increasing area)", y = "Site (by increasing area)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p3)

## BR Heat map of TDN synchrony ##
# Build symmetric df using area-ordered factor levels
brarea_tdn <- bind_rows(
  nresults_br_flagged,
  nresults_br_flagged %>% rename(site_x = site_y, site_y = site_x),
  data.frame(site_x = brsites_by_area, site_y = brsites_by_area, synchrony = 1)
) %>%
  mutate(
    site_x = factor(site_x, levels = brsites_by_area),
    site_y = factor(site_y, levels = brsites_by_area)
  )

p4 <- ggplot(brarea_tdn, aes(x = site_x, y = site_y, fill = synchrony)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(synchrony, 2)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Synchrony"
  ) +
  labs(
    title = "Brush Creek Pairwise TDN Synchrony",
    x = "Site (by increasing area)", y = "Site (by increasing area)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p4)

#### Synchrony plots by watershed area ####
## NM TDN//DOC Synchrony by watershed area ##
# Join area and order sites by increasing watershed area
results_nm_area <- cnresults_nm %>%
  left_join(nmarea_lookup, by = "Site") %>%
  mutate(Site = factor(Site, levels = nmarea_lookup$Site))

p5 <- ggplot(results_nm_area, aes(x = Site, y = synchrony)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "gray40") +
  geom_point(aes(fill = synchrony), shape = 21, size = 4, stroke = 0.5) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Synchrony"
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title = "USF Synchrony between NPOC and TDN",
    subtitle = "Sites ordered by increasing watershed area",
    x = "Site",
    y = "Synchrony"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
print(p5)

## BR TDN//DOC Synchrony by watershed area ##
# Join area and order sites by increasing watershed area
results_br_area <- cnresults_br %>%
  left_join(brarea_lookup, by = "Site") %>%
  mutate(Site = factor(Site, levels = brarea_lookup$Site))

p6 <- ggplot(results_br_area, aes(x = Site, y = synchrony)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", color = "gray40") +
  geom_point(aes(fill = synchrony), shape = 21, size = 4, stroke = 0.5) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Synchrony"
  ) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title = "Brush Creek Synchrony between NPOC and TDN",
    subtitle = "Sites ordered by increasing watershed area",
    x = "Site",
    y = "Synchrony"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
print(p6)




########## Testing the package function!!! ############

#we'll see if we get the same values
#toy datasets are long by site, so we'll use the tsync fncn that reshapes them automatically
#pulling the real package function instead of a hand-copy, so this can't drift out of sync again

correlation_matrix_to_pairs <- function(correlations, data_matrix) {
  sites <- rownames(correlations)
  idx <- which(upper.tri(correlations), arr.ind = TRUE)
  data.frame(
    site_x      = sites[idx[, "row"]],
    site_y      = sites[idx[, "col"]],
    correlation = correlations[idx],
    n_shared    = apply(idx, 1, function(ij) sum(stats::complete.cases(data_matrix[, ij]))),
    row.names   = NULL
  )
}

n_shared_matrix <- function(data_matrix) {
  n <- ncol(data_matrix)
  shared <- matrix(NA_integer_, n, n, dimnames = list(colnames(data_matrix), colnames(data_matrix)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      shared[i, j] <- sum(stats::complete.cases(data_matrix[, c(i, j)]))
    }
  }
  shared
}

tsync <- function(data, conc_col, site_col, time_col, method = c("spearman", "pearson"), min_shared = NULL) {
  method <- match.arg(method)
  df <- data[, c(time_col, site_col, conc_col)]
  df[[conc_col]] <- as.numeric(df[[conc_col]])

  wide <- reshape(df, idvar = time_col, timevar = site_col, direction = "wide")
  wide[[time_col]] <- NULL
  colnames(wide) <- sub(paste0("^", conc_col, "\\."), "", colnames(wide))
  wide <- as.matrix(wide)

  rcorr <- stats::cor(wide, method = method, use = "pairwise.complete.obs")
  diag(rcorr) <- NA

  if (!is.null(min_shared)) {
    rcorr[n_shared_matrix(wide) < min_shared] <- NA_real_
  }

  by_site <- apply(rcorr, 1, stats::median, na.rm = TRUE)

  list(
    by_site = by_site,
    pairwise = correlation_matrix_to_pairs(rcorr, wide)
  )
}

#look at TDN synchrony for brush creek and compare that

site_col = "Site"
time_col = "Date"
conc_col = "TDN..mg.N.L."

sync_values <- tsync(br_toy, conc_col, site_col, time_col, method = "pearson", min_shared = 0)

print('Temporal Synchrony Value pairs:')
print(sync_values$pairwise)

print('Temporal Synchrony Values by site:')
print(sync_values$by_site)

average_sync <- mean(sync_values$by_site, na.rm = TRUE)
print(paste('Average Temporal Synchrony:', average_sync))



########## Recreating an example figure using tsync() ##########

# min_shared is now built into tsync() itself, so the "< 3 shared dates -> NA" flagging
# the original heatmaps use just falls out of the call directly.
nresults_br_pkg_flagged <- tsync(br_toy, "TDN..mg.N.L.", "Site", "Date", method = "pearson", min_shared = 3)$pairwise %>%
  rename(synchrony = correlation)


#### BR Heat map of TDN synchrony (tsync) ####
brarea_tdn_pkg <- bind_rows(
  nresults_br_pkg_flagged,
  nresults_br_pkg_flagged %>% rename(site_x = site_y, site_y = site_x),
  data.frame(site_x = brsites_by_area, site_y = brsites_by_area, synchrony = 1)
) %>%
  mutate(
    site_x = factor(site_x, levels = brsites_by_area),
    site_y = factor(site_y, levels = brsites_by_area)
  )

p4_pkg <- ggplot(brarea_tdn_pkg, aes(x = site_x, y = site_y, fill = synchrony)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(synchrony, 2)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Synchrony"
  ) +
  labs(
    title = "Brush Creek Pairwise TDN Synchrony (tsync package function)",
    x = "Site (by increasing area)", y = "Site (by increasing area)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p4_pkg)


##########################################
#### Monte Carlo sensitivity analysis ####
##########################################
#### (NM) Monte Carlo sensitivity analysis: temporal sample size and synchrony ####
# For three synchrony metrics -- spatial synchrony of NPOC between 2 sites (between-location correlation),
# spatial synchrony of TDN between 2 sites, and cross-variable synchrony between NPOC and TDN at 1 sites
# -- draw n_iter independent synthetic draws at each campaign count in
# camp_grid, and see how the resulting synchrony estimate's median and spread compare to Ref_Sync (the
# "full 5-year record" truth from the static extended dataset).

# Each draw does the same four things:
# 1. pick n random campaign dates across the 5-year synthetic span (shared by both constituents, so the
#    same simulated campaigns can be used for the cross-variable synchrony too)
# 2. build spatially-correlated random noise for those n campaigns (using SSN2's fit, L_space below), for
#    each constituent; the two constituents' noise draws are correlated at rho_cross before the spatial step
# 3. correlate that noise across the n campaigns (red noise), based on how many days apart each pair of
#    them is, using the fitted persistence estimate (rho)
# 4. add SSN2's spatial mean and seasonal cycle, exponentiate back to the natural (non-log) scale
# The loops below iterate over campaignn and n_iter random draws

# Correlation between the two constituents' noise draws at the same location and time, before
# spatial/temporal structure is applied. 0 = no shared noise (the two constituents are linked only through
# shared seasonality). Set this to match how synthetic_extended's cross-constituent correlation was built,
# if that's known; otherwise leave at 0 and treat the Cross panel as a lower bound on how fast
# cross-variable synchrony stabilizes.
rho_cross <- 0


if (length(constituents) != 2) {
  stop("This synchrony Monte Carlo expects exactly two constituents (spatial synchrony for each, plus ",
       "cross-variable synchrony between them); constituents has length ", length(constituents), ".")
}

# n campaigns below which Pearson r is degenerate (|r| = 1 at n = 2); mirrors the < 3 shared dates flagging
# used elsewhere in this script.
sync_grid <- camp_grid[camp_grid >= 3]

mean_offdiag <- function(cm) mean(cm[upper.tri(cm)], na.rm = TRUE)

# Row-wise Pearson correlation between two matrices (row i of A vs row i of B)
row_cor <- function(A, B) {
  Ac <- A - rowMeans(A)
  Bc <- B - rowMeans(B)
  rowSums(Ac * Bc) / sqrt(rowSums(Ac^2) * rowSums(Bc^2))
}

## Reference synchrony: the static extended dataset's spatial and cross-variable correlations ##
# predID x CampaignNum matrix of Conc, one per constituent
wide_extended <- lapply(constituents, function(cc) {
  nm_toy_extended %>%
    filter(Constituent == cc) %>%
    select(predID, CampaignNum, Conc) %>%
    pivot_wider(names_from = CampaignNum, values_from = Conc) %>%
    arrange(predID)
})
names(wide_extended) <- constituents

# Keep only locations shared by both constituents, in the same order
common_ids    <- Reduce(intersect, lapply(wide_extended, function(w) w$predID))
wide_extended <- lapply(wide_extended, function(w) {
  as.matrix(w[match(common_ids, w$predID), setdiff(names(w), "predID")])
})

# Reference spatial synchrony: the static extended dataset's between-location correlation, across all 60
# synthetic campaigns, averaged across all location pairs. One value per constituent.
ref_sync_spatial <- sapply(constituents, function(cc) mean_offdiag(cor(t(wide_extended[[cc]]))))
names(ref_sync_spatial) <- constituents

# Reference cross-variable synchrony: the static extended dataset's constituents[1]-constituents[2]
# correlation across all 60 synthetic campaigns at each location, averaged across all m locations.
ref_sync_cross <- mean(sapply(seq_along(common_ids), function(j) {
  cor(wide_extended[[constituents[1]]][j, ], wide_extended[[constituents[2]]][j, ])
}), na.rm = TRUE)

ref_sync     <- c(ref_sync_spatial, Cross = ref_sync_cross)
metric_names <- c(constituents, "Cross")

## Field setup, shared m, validated rho ##
fs <- nm_field_setups[constituents]   # this constituent pair's two field setups

for (cc in constituents) {
  rho <- fs[[cc]]$temporal_rho
  if (rho < 0 || rho >= 1) {
    stop(cc, " temporal_rho (", rho, ") is outside [0, 1) -- this exponential-decay ",
         "temporal correlation can't represent a negative value; see 02_build_synthetic_data.R's Part D.")
  }
}

m <- fs[[constituents[1]]]$setup$m
if (fs[[constituents[2]]]$setup$m != m) {
  stop("field_setups for ", constituents[1], " and ", constituents[2],
       " have different numbers of predpts locations -- can't compute cross-variable synchrony.")
}

# One synthetic field (m x n, natural scale) for one constituent, given a standard-normal draw z (m x n)
simulate_field <- function(cc, z, t_sample, doy_sample) {
  setup   <- fs[[cc]]$setup
  beta    <- fs[[cc]]$beta
  rho     <- fs[[cc]]$temporal_rho
  ref_doy <- fs[[cc]]$ref_doy
  n       <- length(t_sample)
  
  # spatially-correlated noise, then correlated across campaigns (red noise)
  spatial_noise <- t(setup$L_space) %*% z
  elapsed_days  <- abs(outer(t_sample, t_sample, "-")) / time_unit_days
  time_chol     <- chol(rho ^ elapsed_days + diag(1e-8, n))
  noise         <- spatial_noise %*% t(time_chol)
  
  # seasonal shift, relative to ref_doy
  ref_shift <- beta[3] * sin(2 * pi * ref_doy / 365) + beta[4] * cos(2 * pi * ref_doy / 365)
  shift     <- beta[3] * sin(2 * pi * doy_sample / 365) + beta[4] * cos(2 * pi * doy_sample / 365) - ref_shift
  shift_mat <- matrix(shift, nrow = m, ncol = n, byrow = TRUE)
  
  exp(setup$mu_cond + shift_mat + noise)  # m x n, natural scale
}

## The Monte Carlo loop ##
sync_mc_results_list <- list()

for (n in sync_grid) {
  
  draws <- matrix(NA_real_, nrow = n_iter, ncol = length(metric_names),
                  dimnames = list(NULL, metric_names))
  
  for (i in seq_len(n_iter)) {
    
    # Step 1: draw n random campaign dates across the 5-year synthetic span, shared by both constituents
    t_sample   <- sort(runif(n, 0, total_days))
    doy_sample <- t_sample %% 365.25
    
    # Step 2: standard-normal draws, correlated across the two constituents at rho_cross
    z1 <- matrix(rnorm(m * n), nrow = m, ncol = n)
    z2 <- rho_cross * z1 + sqrt(1 - rho_cross^2) * matrix(rnorm(m * n), nrow = m, ncol = n)
    
    # Steps 3-4: spatial + temporal structure, seasonal shift, back-transform
    field1 <- simulate_field(constituents[1], z1, t_sample, doy_sample)
    field2 <- simulate_field(constituents[2], z2, t_sample, doy_sample)
    
    # Synchrony metrics for this draw
    draws[i, constituents[1]] <- mean_offdiag(cor(t(field1)))
    draws[i, constituents[2]] <- mean_offdiag(cor(t(field2)))
    draws[i, "Cross"]         <- mean(row_cor(field1, field2), na.rm = TRUE)
  }
  
  sync_mc_results_list[[paste(n)]] <- as_tibble(draws) %>%
    mutate(N = n) %>%
    pivot_longer(all_of(metric_names), names_to = "Metric", values_to = "Sync") %>%
    mutate(Ref_Sync = ref_sync[Metric])
  
  message("Finished N = ", n)
}

## Combine, summarize, and plot
sync_mc_results <- bind_rows(sync_mc_results_list)

sync_mc_summary <- sync_mc_results %>%
  group_by(Metric, N) %>%
  summarize(Median_Sync = median(Sync), P05_Sync = quantile(Sync, 0.05), P95_Sync = quantile(Sync, 0.95),
            Ref_Sync = first(Ref_Sync), .groups = "drop") %>%
  mutate(Bias = Median_Sync - Ref_Sync)

print(sync_mc_summary, n = Inf)
write_csv(sync_mc_summary, file.path(data_out_dir, "synchrony_MC_summary.csv"))

# One row per metric, this analysis's own Bias at the N closest to actual_n_campaigns (sync_grid may not
# contain actual_n_campaigns exactly, since counts below 3 are dropped)
bias_annot <- sync_mc_summary %>%
  group_by(Metric) %>%
  slice_min(abs(N - actual_n_campaigns), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = sprintf("At N=%d: %+.2f bias", N, Bias))

p_sync <- ggplot(sync_mc_summary, aes(x = N, y = Median_Sync)) +
  geom_ribbon(aes(ymin = P05_Sync, ymax = P95_Sync), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_Sync), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = actual_n_campaigns, color = "red", linetype = "dashed") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  labs(x = "Sample size (N campaigns)", y = "Synchrony (median, 5th-95th pct band)",
       title = "Temporal sample-size sensitivity of NM synchrony",
       subtitle = "Blue = reference synchrony (static extended dataset, all 60 synthetic months); red = actual sample size, annotated with bias at the nearest tested N. Each Monte Carlo iteration redraws a fresh spatiotemporal field, campaign times drawn at random across a 5-year span.") +
  theme_bw()

print(p_sync)
ggsave(file.path(plot_dir, "synchrony_MC_plot_nm.png"), p_sync, width = 10, height = 4, dpi = 150)

#### (BR) Monte Carlo sensitivity analysis: temporal sample size and synchrony ####
# For three synchrony metrics -- spatial synchrony of NPOC between 2 sites (between-location correlation),
# spatial synchrony of TDN between 2 sites, and cross-variable synchrony between NPOC and TDN at 1 sites
# -- draw n_iter independent synthetic draws at each campaign count in
# camp_grid, and see how the resulting synchrony estimate's median and spread compare to Ref_Sync (the
# "full 5-year record" truth from the static extended dataset).

# Each draw does the same four things:
# 1. pick n random campaign dates across the 5-year synthetic span (shared by both constituents, so the
#    same simulated campaigns can be used for the cross-variable synchrony too)
# 2. build spatially-correlated random noise for those n campaigns (using SSN2's fit, L_space below), for
#    each constituent; the two constituents' noise draws are correlated at rho_cross before the spatial step
# 3. correlate that noise across the n campaigns (red noise), based on how many days apart each pair of
#    them is, using the fitted persistence estimate (rho)
# 4. add SSN2's spatial mean and seasonal cycle, exponentiate back to the natural (non-log) scale
# The loops below iterate over campaignn and n_iter random draws

# Correlation between the two constituents' noise draws at the same location and time, before
# spatial/temporal structure is applied. 0 = no shared noise (the two constituents are linked only through
# shared seasonality). Set this to match how synthetic_extended's cross-constituent correlation was built,
# if that's known; otherwise leave at 0 and treat the Cross panel as a lower bound on how fast
# cross-variable synchrony stabilizes.
rho_cross <- 0


if (length(constituents) != 2) {
  stop("This synchrony Monte Carlo expects exactly two constituents (spatial synchrony for each, plus ",
       "cross-variable synchrony between them); constituents has length ", length(constituents), ".")
}

# n campaigns below which Pearson r is degenerate (|r| = 1 at n = 2); mirrors the < 3 shared dates flagging
# used elsewhere in this script.
sync_grid <- camp_grid[camp_grid >= 3]

mean_offdiag <- function(cm) mean(cm[upper.tri(cm)], na.rm = TRUE)

# Row-wise Pearson correlation between two matrices (row i of A vs row i of B)
row_cor <- function(A, B) {
  Ac <- A - rowMeans(A)
  Bc <- B - rowMeans(B)
  rowSums(Ac * Bc) / sqrt(rowSums(Ac^2) * rowSums(Bc^2))
}

## Reference synchrony: the static extended dataset's spatial and cross-variable correlations ##
# predID x CampaignNum matrix of Conc, one per constituent
wide_extended <- lapply(constituents, function(cc) {
  br_toy_extended %>%
    filter(Constituent == cc) %>%
    select(predID, CampaignNum, Conc) %>%
    pivot_wider(names_from = CampaignNum, values_from = Conc) %>%
    arrange(predID)
})
names(wide_extended) <- constituents

# Keep only locations shared by both constituents, in the same order
common_ids    <- Reduce(intersect, lapply(wide_extended, function(w) w$predID))
wide_extended <- lapply(wide_extended, function(w) {
  as.matrix(w[match(common_ids, w$predID), setdiff(names(w), "predID")])
})

# Reference spatial synchrony: the static extended dataset's between-location correlation, across all 60
# synthetic campaigns, averaged across all location pairs. One value per constituent.
ref_sync_spatial <- sapply(constituents, function(cc) mean_offdiag(cor(t(wide_extended[[cc]]))))
names(ref_sync_spatial) <- constituents

# Reference cross-variable synchrony: the static extended dataset's constituents[1]-constituents[2]
# correlation across all 60 synthetic campaigns at each location, averaged across all m locations.
ref_sync_cross <- mean(sapply(seq_along(common_ids), function(j) {
  cor(wide_extended[[constituents[1]]][j, ], wide_extended[[constituents[2]]][j, ])
}), na.rm = TRUE)

ref_sync     <- c(ref_sync_spatial, Cross = ref_sync_cross)
metric_names <- c(constituents, "Cross")

## Field setup, shared m, validated rho ##
fs <- br_field_setups[constituents]   # this constituent pair's two field setups

for (cc in constituents) {
  rho <- fs[[cc]]$temporal_rho
  if (rho < 0 || rho >= 1) {
    stop(cc, " temporal_rho (", rho, ") is outside [0, 1) -- this exponential-decay ",
         "temporal correlation can't represent a negative value; see 02_build_synthetic_data.R's Part D.")
  }
}

m <- fs[[constituents[1]]]$setup$m
if (fs[[constituents[2]]]$setup$m != m) {
  stop("field_setups for ", constituents[1], " and ", constituents[2],
       " have different numbers of predpts locations -- can't compute cross-variable synchrony.")
}

# One synthetic field (m x n, natural scale) for one constituent, given a standard-normal draw z (m x n)
simulate_field <- function(cc, z, t_sample, doy_sample) {
  setup   <- fs[[cc]]$setup
  beta    <- fs[[cc]]$beta
  rho     <- fs[[cc]]$temporal_rho
  ref_doy <- fs[[cc]]$ref_doy
  n       <- length(t_sample)
  
  # spatially-correlated noise, then correlated across campaigns (red noise)
  spatial_noise <- t(setup$L_space) %*% z
  elapsed_days  <- abs(outer(t_sample, t_sample, "-")) / time_unit_days
  time_chol     <- chol(rho ^ elapsed_days + diag(1e-8, n))
  noise         <- spatial_noise %*% t(time_chol)
  
  # seasonal shift, relative to ref_doy
  ref_shift <- beta[3] * sin(2 * pi * ref_doy / 365) + beta[4] * cos(2 * pi * ref_doy / 365)
  shift     <- beta[3] * sin(2 * pi * doy_sample / 365) + beta[4] * cos(2 * pi * doy_sample / 365) - ref_shift
  shift_mat <- matrix(shift, nrow = m, ncol = n, byrow = TRUE)
  
  exp(setup$mu_cond + shift_mat + noise)  # m x n, natural scale
}

## The Monte Carlo loop ##
sync_mc_results_list <- list()

for (n in sync_grid) {
  
  draws <- matrix(NA_real_, nrow = n_iter, ncol = length(metric_names),
                  dimnames = list(NULL, metric_names))
  
  for (i in seq_len(n_iter)) {
    
    # Step 1: draw n random campaign dates across the 5-year synthetic span, shared by both constituents
    t_sample   <- sort(runif(n, 0, total_days))
    doy_sample <- t_sample %% 365.25
    
    # Step 2: standard-normal draws, correlated across the two constituents at rho_cross
    z1 <- matrix(rnorm(m * n), nrow = m, ncol = n)
    z2 <- rho_cross * z1 + sqrt(1 - rho_cross^2) * matrix(rnorm(m * n), nrow = m, ncol = n)
    
    # Steps 3-4: spatial + temporal structure, seasonal shift, back-transform
    field1 <- simulate_field(constituents[1], z1, t_sample, doy_sample)
    field2 <- simulate_field(constituents[2], z2, t_sample, doy_sample)
    
    # Synchrony metrics for this draw
    draws[i, constituents[1]] <- mean_offdiag(cor(t(field1)))
    draws[i, constituents[2]] <- mean_offdiag(cor(t(field2)))
    draws[i, "Cross"]         <- mean(row_cor(field1, field2), na.rm = TRUE)
  }
  
  sync_mc_results_list[[paste(n)]] <- as_tibble(draws) %>%
    mutate(N = n) %>%
    pivot_longer(all_of(metric_names), names_to = "Metric", values_to = "Sync") %>%
    mutate(Ref_Sync = ref_sync[Metric])
  
  message("Finished N = ", n)
}

## Combine, summarize, and plot
sync_mc_results <- bind_rows(sync_mc_results_list)

sync_mc_summary <- sync_mc_results %>%
  group_by(Metric, N) %>%
  summarize(Median_Sync = median(Sync), P05_Sync = quantile(Sync, 0.05), P95_Sync = quantile(Sync, 0.95),
            Ref_Sync = first(Ref_Sync), .groups = "drop") %>%
  mutate(Bias = Median_Sync - Ref_Sync)

print(sync_mc_summary, n = Inf)
write_csv(sync_mc_summary, file.path(data_out_dir, "synchrony_MC_summary_br.csv"))

# One row per metric, this analysis's own Bias at the N closest to actual_n_campaigns (sync_grid may not
# contain actual_n_campaigns exactly, since counts below 3 are dropped)
bias_annot <- sync_mc_summary %>%
  group_by(Metric) %>%
  slice_min(abs(N - actual_n_campaigns), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = sprintf("At N=%d: %+.2f bias", N, Bias))

p_sync <- ggplot(sync_mc_summary, aes(x = N, y = Median_Sync)) +
  geom_ribbon(aes(ymin = P05_Sync, ymax = P95_Sync), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_Sync), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = actual_n_campaigns, color = "red", linetype = "dashed") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  labs(x = "Sample size (N campaigns)", y = "Synchrony (median, 5th-95th pct band)",
       title = "Temporal sample-size sensitivity of BR synchrony",
       subtitle = "Blue = reference synchrony (static extended dataset, all 60 synthetic months); red = actual sample size, annotated with bias at the nearest tested N. Each Monte Carlo iteration redraws a fresh spatiotemporal field, campaign times drawn at random across a 5-year span.") +
  theme_bw()

print(p_sync)
ggsave(file.path(plot_dir, "synchrony_MC_plot_br.png"), p_sync, width = 10, height = 4, dpi = 150)
