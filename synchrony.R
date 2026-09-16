#### Read me ####
# This script is used to calculate synchrony between sites within a watershed with a toy dataset for QuEST's planned commentary manuscript on the use of spatiotemporal metrics for assessing the spatiotemporal variance of stream chemistry across stream networks.
# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Eva Tipps, 2026-09-04 (with help building functions from Claude version 1.24012.9 (03c61d) 2026-09-04T09:58:04.000Z... reviewed and edited by E. Tipps)
# Last update (Person, Date): Eva Tipps, 2026-09-16

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
nm_field_setups     <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
nm_toy_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)

br_toy           <- read_csv(file.path(data_out_dir, "br_clean.csv"), show_col_types = FALSE)
br_field_setups     <- readRDS(file.path(data_out_dir, "br_field_setups.rds"))
br_toy_extended <- read_csv(file.path(data_out_dir, "br_synthetic_extended.csv"), show_col_types = FALSE)


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

#### Synchrony between each site for NPOC ####
# NM #
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

# BR #
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



#### Synchrony between each site for TDN ####
# NM #
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

# BR #
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


#### Synchrony between NPOC and TDN at each site ####
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

#### NM Heat map of DOC synchrony ####
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

#### BR Heat map of DOC synchrony ####
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


#### NM Heat map of TDN synchrony ####

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

#### BR Heat map of TDN synchrony ####
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


#### NM TDN//DOC Synchrony by watershed area ####
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

#### BR TDN//DOC Synchrony by watershed area ####
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

