#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-08-23 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-08
# Bre Rivera Waterman, 2026-09-09 pulling in package/helper function and comparing to previous calculations
# Bre Rivera Waterman,09/29/2026 for conceptual fig values
# Alex Webster, 2026-10-06 to add br analysis and remove sensitivity analysis
# Alex Webster, 2026-10-08 to add extended synthetic data to make CVt n (3 per watershed) closer to CVs n (~20) for comparison


# Requires: 02_build_synthetic_data.R must be run first to produce data/[dataset]_clean.csv and data/[dataset]_synthetic_extended.csv. If these are available, no need to rerun 02_build_synthetic_data.R

# This script demonstrates how to calculate CVt, the temporal coefficient of variation (across sampling campaigns, one value per sampling site) on three datasets: NM, BR, and Gu et al. 2021. 
# Complementary sensitivity analyses can be found in CVt_sensitivity.R

# Outputs:
#   data/CVt_observed.csv:      real CVt per catchment x site x constituent
#   data/CVt_observed_mean.csv: real CVt averaged across catchment x sites, per constituent
#   data/CVt_observed_extended.csv:      CVt using extended synthetic data per catchment x site x constituent
#   data/CVt_observed_extended_mean.csv: CVt using extended synthetic data, averaged across catchment x sites, per constituent

#### Packages ####
library(tidyverse)
source("cv_fxhelper.R")
library(SSN2)

#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

# chem data from 3 consecutive real campaigns in nm and br
clean_nm            <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
clean_nm$Catchment = "nm"
clean_br            <- read_csv(file.path(data_out_dir, "br_clean.csv"), show_col_types = FALSE)
clean_br$Catchment = "br"

clean <- bind_rows(clean_nm, clean_br)

actual_n_campaigns <- n_distinct(clean$CampaignID)  # real number of campaigns conducted across catchments
actual_n_campaigns_nm <- n_distinct(clean_nm$CampaignID)  # real number of nm campaigns conducted
actual_n_campaigns_br <- n_distinct(clean_br$CampaignID)  # real number of br campaigns conducted

# chem data from synthetically extended data to 20 campaigns
extended_nm            <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)
extended_nm$Catchment = "nm"
#extended_nm = extended_nm %>% filter(CampaignNum %in% c(1:20)) %>% filter(predID %in% c())
extended_br            <- read_csv(file.path(data_out_dir, "br_synthetic_extended.csv"), show_col_types = FALSE)
extended_br$Catchment = "br"
#extended_br = extended_br %>% filter(CampaignNum %in% c(1:20))

extended <- bind_rows(extended_nm, extended_br)

extended_n_campaigns <- n_distinct(extended$CampaignNum)  # number of campaigns across catchments
extended_n_campaigns_nm <- n_distinct(extended_nm$CampaignNum)  # number of nm campaigns 
extended_n_campaigns_br <- n_distinct(extended_br$CampaignNum)  # number of br campaigns 

#### PART A.1 -- Calculate CVt of real toy dataset ####
# Per-site CV across that site's real campaigns

## For each constituent in `constituents`: per-site CV across that site's real campaigns (sites with fewer than 3 real campaigns are dropped, since sd() of 1-2 points is too noisy to call a temporal CV).
CVt_observed_list      <- list()
CVt_observed_mean_list <- list()

for (cc in constituents) {
  
  per_site_cvt <- clean %>%
    filter(!is.na(.data[[cc]])) %>%
    group_by(Catchment, Site) %>%
    summarize(n_sites = n(), CVt = sd(.data[[cc]]) / mean(.data[[cc]]), .groups = "drop") %>%
    filter(n_sites >= 3)

  print(per_site_cvt)
  
  CVt_observed_list[[cc]] <- per_site_cvt %>%
    mutate(Constituent = cc, .before = 1)
  
  CVt_observed_mean_list[[cc]] <- per_site_cvt %>%
    group_by(Catchment) %>%
    summarize(n_sites_used = n(),
              CVt_observed = mean(CVt),
              .groups = "drop") %>%
    mutate(Constituent = cc, .before = 1)
}


## ---- Combine and save ----
CVt_observed <- bind_rows(CVt_observed_list)
write_csv(CVt_observed, file.path(data_out_dir, "CVt_observed.csv"))

CVt_observed_mean <- bind_rows(CVt_observed_mean_list)
write_csv(CVt_observed_mean, file.path(data_out_dir, "CVt_observed_mean.csv"))
print(CVt_observed_mean)


#### PART A.2 -- Calculate CVt of real toy dataset using package!! ####
# CVt is calculated for one site across campaigns.

CVt_result <- temporal_cv(data = clean,
                         concentration = constituents,
                         site = "Site",
                         watershed = "Catchment",
                         digits = 2)

# One row per site × constituent
CVt_by_site.2 <- 
  CVt_result$by_site %>%
  transmute(Catchment = Catchment,
            Constituent = constituent,
            Site,
            n_sites = n_used,
            CVt = temporal_cv)

print(CVt_by_site.2)

write_csv(CVt_by_site.2, file.path(data_out_dir, "CVt_by_site.2.csv"))

# Mean and SD across campaign-level CVs for each constituent
CVt_observed.2 <- CVt_result$watershed_summary %>%
  transmute(Catchment = Catchment,
            Constituent = constituent,
            n_sites = n_sites_used,
            Mean_CVt = mean_temporal_cv,
            SD_CVt = sd_temporal_cv )

print(CVt_observed.2)

write_csv(CVt_observed.2, file.path(data_out_dir, "CVt_observed.2.csv"))

#### PART A.3 -- Compare original and package results ####
site_comparison <- CVt_observed %>%
  mutate(CVt_original_rounded = round(CVt, 2)) %>%
  full_join(CVt_by_site.2 %>%
              rename(
                #n_sites_package = n_sites,
                     CVt_package = CVt ),
            by = c("Catchment","Constituent", "Site")) %>%
  #rename(n_sites_original = n_sites) %>%
  mutate(#n_sites_difference = n_sites_package - n_sites_original,
         CVt_difference = CVt_package - CVt_original_rounded) %>%
  arrange(Constituent, Site)

print(site_comparison)


#watershed summaries
original_by_catchment <- CVt_observed %>%
  group_by(Constituent, Catchment) %>%
  summarize(n_sites_original = n(),
            Mean_CVt_original = mean(CVt),
            SD_CVt_original = if (n() >= 2) sd(CVt) else NA_real_,
            .groups = "drop") %>%
  mutate(Mean_CVt_original_rounded = round(Mean_CVt_original, 2),
         SD_CVt_original_rounded = round(SD_CVt_original, 2))

package_by_catchment <- CVt_observed.2 %>%
  rename(n_sites_package = n_sites,
         Mean_CVt_package = Mean_CVt,
         SD_CVt_package = SD_CVt)

summary_comparison <- original_by_catchment %>%
  full_join(package_by_catchment, by = c("Constituent", "Catchment")) %>%
  mutate(Mean_CVt_difference = Mean_CVt_package - Mean_CVt_original_rounded,
         SD_CVt_difference = SD_CVt_package - SD_CVt_original_rounded) %>%
  arrange(Constituent, Catchment)

print(summary_comparison, n = Inf)


#visual comparison
CVt_plot_data <- 
  bind_rows(CVt_observed %>%
              transmute(Catchment, Constituent, Site, CVt,
                        Approach = "Original QuEST calculation" ),
            CVt_by_site.2 %>%
              transmute(Catchment, Constituent,Site, CVt, Approach = "CV helper")) %>%
  filter(!is.na(CVt)) %>%
  mutate(Approach = factor(Approach,
                           levels = c("Original QuEST calculation", "CV helper") ) )


p_cvt_comparison <- 
  ggplot(CVt_plot_data, aes(x = Approach, y = CVt, fill = Approach)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45,
    color = NA ) +
  geom_jitter(
    aes(color = Catchment),
    width = 0.05,
    height = 0,
    size = 1.8,
    alpha = 0.75 ) +
  facet_wrap(~ Constituent,
             scales = "free_y",
             ncol = 2) +
  scale_fill_manual(
    values = c(
      "Original QuEST calculation" = "#4C78A8",
      "CV helper" = "#F58518" )) +
  labs(title = "Temporal CV by calculation approach",
       subtitle = "Each point is one sampling location",
       x = NULL,
       y = "Temporal coefficient of variation (CVt)",
       fill = "Approach") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

p_cvt_comparison

#### PART A.4 -- calculate using extended synthetic data ####

### the purpose of this section is to add extended synthetic data to make CVt n (3 per watershed) closer to CVs n (~20) for comparison.
# max CV is a function of n (CVmax = √n), so CVt and CVs having different n impacts the comparison of their relative values and thus if you're able to accurately interpret whether spatial or temporal variance dominates in a given dataset. We are going to side-step this issue by matching their samples sizes, and test sensitivities to this further in CVt_sensitivity.R.

### First, compare the real and synthetic data equivalent in time and space. 
# This is to assess how to use the synthetic data to increase # of campaigns while keeping it as similar as possible to the real data results. 
# In the real data, CV nm > CV br for both DOC and TDN and both CVt and CVs.
# Need to decide whether I'm using synthetic data for just CVt or both CVt and CVs, or some combination of real and synthetic data. 

circ_dist <- function(a, b) pmin(abs(a - b), 365.25 - abs(a - b))  # shared helper, not catchment-specific

## nm - spatial match
ssn_object_nm <- readRDS("data/ssn_objects.rds")$nm
dist_list_nm  <- ssn_get_stream_distmat(ssn_object_nm, name = "predpts")
# total routed distance, predpts (rows) x obs (cols):
total_dist_nm <- t(dist_list_nm$dist.net1.a) + dist_list_nm$dist.net1.b
# pid <-> Site/predID lookups
obs_geom_nm  <- ssn_get_netgeom(ssn_object_nm$obs, reformat = TRUE)
obs_geom_nm$Site <- ssn_object_nm$obs$Site
pred_geom_nm <- ssn_get_netgeom(ssn_object_nm$preds$predpts, reformat = TRUE)
pred_geom_nm$predID <- ssn_object_nm$preds$predpts$predID
site_lookup_nm <- setNames(obs_geom_nm$Site, as.character(obs_geom_nm$pid))     # pid -> Site
pred_lookup_nm <- setNames(pred_geom_nm$predID, as.character(pred_geom_nm$pid)) # pid -> predID
# nearest predpt per real site
a_nm <- t(dist_list_nm$dist.net1.a)   # pred (rows) x obs (cols): >0 iff pred is upstream of obs
b_nm <- dist_list_nm$dist.net1.b      # pred (rows) x obs (cols): >0 iff pred is downstream of obs
upstream_dist_nm <- a_nm
upstream_dist_nm[!(a_nm > 0 & b_nm == 0)] <- NA   # keep ONLY true upstream, flow-connected pairs- drops flow-unconnected pairs (both >0)
nearest_row_up_nm  <- apply(upstream_dist_nm, 2, function(col) if (all(is.na(col))) NA_integer_ else which.min(col))
nearest_dist_up_nm <- apply(upstream_dist_nm, 2, function(col) if (all(is.na(col))) NA_real_    else min(col, na.rm = TRUE))
site_match_nm_upstream <- tibble(
  Site               = site_lookup_nm[colnames(upstream_dist_nm)],
  nearest_predID     = ifelse(is.na(nearest_row_up_nm), NA_character_,
                              pred_lookup_nm[rownames(upstream_dist_nm)[nearest_row_up_nm]]),
  network_distance_m = nearest_dist_up_nm
)
sum(is.na(site_match_nm_upstream$nearest_predID))
## nm - temporal match
campaign_doy_nm <- clean_nm %>% group_by(CampaignID) %>% summarize(DOY = mean(DOY), .groups = "drop")
synth_doy_nm    <- extended_nm %>% distinct(CampaignNum, DOY)
# circular DOY distance, synthetic months (rows) x real campaigns (cols)
doy_dist_nm <- outer(synth_doy_nm$DOY, campaign_doy_nm$DOY, FUN = circ_dist)
rownames(doy_dist_nm) <- synth_doy_nm$CampaignNum
colnames(doy_dist_nm) <- campaign_doy_nm$CampaignID
tie_tol_days     <- 1
campaign_nums_nm <- as.integer(rownames(doy_dist_nm))
pick_nearest_campaign <- function(col, campaign_nums, tol) {
  best_gap       <- min(col)
  idx_candidates <- which(col <= best_gap + tol)                         # all near-tied rows
  winner         <- idx_candidates[which.min(campaign_nums[idx_candidates])]  # smallest CampaignNum among them
  unname(c(campaign_nums[winner], col[winner]))   # strip names -- col[winner] otherwise carries its
  # own rowname and c() would glue it onto ours
}
picked_nm <- apply(doy_dist_nm, 2, pick_nearest_campaign, campaign_nums = campaign_nums_nm, tol = tie_tol_days)
rownames(picked_nm) <- c("CampaignNum", "day_gap")  # name the rows ourselves now that they're unambiguous
campaign_match_nm <- tibble(
  CampaignID          = colnames(doy_dist_nm),
  nearest_CampaignNum = as.integer(picked_nm["CampaignNum", ]),
  day_gap             = picked_nm["day_gap", ]
)
## nm - build the comparison and summarize
unmatched_sites_nm <- site_match_nm_upstream %>% filter(is.na(nearest_predID)) %>% pull(Site)
if (length(unmatched_sites_nm) > 0) {
  cat("NM: dropping", length(unmatched_sites_nm), "site(s) with no flow-connected predpt upstream:",
      paste(unmatched_sites_nm, collapse = ", "), "\n")
}
comparison_nm <- clean_nm %>%
  select(Site, CampaignID, DOY, all_of(constituents)) %>%
  pivot_longer(all_of(constituents), names_to = "Constituent", values_to = "Conc_real") %>%
  filter(!is.na(Conc_real)) %>%
  inner_join(site_match_nm_upstream %>% filter(!is.na(nearest_predID)), by = "Site") %>%
  inner_join(campaign_match_nm, by = "CampaignID") %>%
  inner_join(extended_nm %>% select(predID, CampaignNum, Constituent, Conc_synth = Conc),
             by = c("nearest_predID" = "predID", "nearest_CampaignNum" = "CampaignNum", "Constituent")) %>%
  mutate(Log_Diff = log(Conc_real) - log(Conc_synth),
         Ratio    = Conc_real / Conc_synth,
         Pct_Diff = (Conc_real - Conc_synth) / Conc_synth * 100)
comparison_summary_nm <- comparison_nm %>%
  group_by(Constituent) %>%
  summarize(n                       = n(),
            Mean_Log_Diff           = mean(Log_Diff),
            SD_Log_Diff             = sd(Log_Diff),
            Median_Pct_Diff         = median(Pct_Diff),
            Mean_network_distance_m = mean(network_distance_m),
            Max_network_distance_m  = max(network_distance_m),
            Mean_day_gap            = mean(day_gap),
            Max_day_gap             = max(day_gap),
            .groups = "drop")
print(comparison_summary_nm)
print(comparison_nm %>% filter(Constituent == "NPOC..mg.C.L.") %>% select(CampaignID, Site, nearest_predID, Conc_real, Conc_synth, Log_Diff), n=50)
print(comparison_nm %>% filter(Constituent == "TDN..mg.N.L.") %>% select(CampaignID, Site, nearest_predID, Conc_real, Conc_synth, Log_Diff), n=50)
# Visual check: real vs. synthetic concentration, one point per site/campaign
ggplot(comparison_nm, aes(x = Conc_real, y = Conc_synth, color = factor(CampaignID))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  #scale_x_log10() +
  #scale_y_log10() +
  facet_wrap(~ Constituent, scales = "free") +
  labs(x = "Real concentration", y = "Synthetic concentration",
       color = "Campaign", title = "NM: real vs. synthetic concentration by site") +
  theme_bw()

## br - spatial match
ssn_object_br <- readRDS("data/ssn_objects_br.rds")$br
dist_list_br  <- ssn_get_stream_distmat(ssn_object_br, name = "predpts")
# total routed distance, predpts (rows) x obs (cols):
total_dist_br <- t(dist_list_br$dist.net1.a) + dist_list_br$dist.net1.b
# pid <-> Site/predID lookups
obs_geom_br  <- ssn_get_netgeom(ssn_object_br$obs, reformat = TRUE)
obs_geom_br$Site <- ssn_object_br$obs$Site
pred_geom_br <- ssn_get_netgeom(ssn_object_br$preds$predpts, reformat = TRUE)
pred_geom_br$predID <- ssn_object_br$preds$predpts$predID
site_lookup_br <- setNames(obs_geom_br$Site, as.character(obs_geom_br$pid))     # pid -> Site
pred_lookup_br <- setNames(pred_geom_br$predID, as.character(pred_geom_br$pid)) # pid -> predID
# nearest predpt per real site
a_br <- t(dist_list_br$dist.net1.a)   # pred (rows) x obs (cols): >0 iff pred is upstream of obs
b_br <- dist_list_br$dist.net1.b      # pred (rows) x obs (cols): >0 iff pred is downstream of obs
upstream_dist_br <- a_br
upstream_dist_br[!(a_br > 0 & b_br == 0)] <- NA   # keep ONLY true upstream, flow-connected pairs- drops flow-unconnected pairs (both >0)
nearest_row_up_br  <- apply(upstream_dist_br, 2, function(col) if (all(is.na(col))) NA_integer_ else which.min(col))
nearest_dist_up_br <- apply(upstream_dist_br, 2, function(col) if (all(is.na(col))) NA_real_    else min(col, na.rm = TRUE))
site_match_br_upstream <- tibble(
  Site               = site_lookup_br[colnames(upstream_dist_br)],
  nearest_predID     = ifelse(is.na(nearest_row_up_br), NA_character_,
                              pred_lookup_br[rownames(upstream_dist_br)[nearest_row_up_br]]),
  network_distance_m = nearest_dist_up_br
)
sum(is.na(site_match_br_upstream$nearest_predID))
## br - temporal match
campaign_doy_br <- clean_br %>% group_by(CampaignID) %>% summarize(DOY = mean(DOY), .groups = "drop")
synth_doy_br    <- extended_br %>% distinct(CampaignNum, DOY)
# circular DOY distance, synthetic months (rows) x real campaigns (cols)
doy_dist_br <- outer(synth_doy_br$DOY, campaign_doy_br$DOY, FUN = circ_dist)
rownames(doy_dist_br) <- synth_doy_br$CampaignNum
colnames(doy_dist_br) <- campaign_doy_br$CampaignID
tie_tol_days     <- 1
campaign_nums_br <- as.integer(rownames(doy_dist_br))
pick_nearest_campaign <- function(col, campaign_nums, tol) {
  best_gap       <- min(col)
  idx_candidates <- which(col <= best_gap + tol)                         # all near-tied rows
  winner         <- idx_candidates[which.min(campaign_nums[idx_candidates])]  # smallest CampaignNum among them
  unname(c(campaign_nums[winner], col[winner]))   # strip names -- col[winner] otherwise carries its
  # own rowname and c() would glue it onto ours
}
picked_br <- apply(doy_dist_br, 2, pick_nearest_campaign, campaign_nums = campaign_nums_br, tol = tie_tol_days)
rownames(picked_br) <- c("CampaignNum", "day_gap")  # name the rows ourselves now that they're unambiguous
campaign_match_br <- tibble(
  CampaignID          = colnames(doy_dist_br),
  nearest_CampaignNum = as.integer(picked_br["CampaignNum", ]),
  day_gap             = picked_br["day_gap", ]
)
## br - build the comparison and summarize
ubratched_sites_br <- site_match_br_upstream %>% filter(is.na(nearest_predID)) %>% pull(Site)
if (length(ubratched_sites_br) > 0) {
  cat("br: dropping", length(ubratched_sites_br), "site(s) with no flow-connected predpt upstream:",
      paste(ubratched_sites_br, collapse = ", "), "\n")
}
comparison_br <- clean_br %>%
  select(Site, CampaignID, DOY, all_of(constituents)) %>%
  pivot_longer(all_of(constituents), names_to = "Constituent", values_to = "Conc_real") %>%
  filter(!is.na(Conc_real)) %>%
  inner_join(site_match_br_upstream %>% filter(!is.na(nearest_predID)), by = "Site") %>%
  inner_join(campaign_match_br, by = "CampaignID") %>%
  inner_join(extended_br %>% select(predID, CampaignNum, Constituent, Conc_synth = Conc),
             by = c("nearest_predID" = "predID", "nearest_CampaignNum" = "CampaignNum", "Constituent")) %>%
  mutate(Log_Diff = log(Conc_real) - log(Conc_synth),
         Ratio    = Conc_real / Conc_synth,
         Pct_Diff = (Conc_real - Conc_synth) / Conc_synth * 100)
comparison_summary_br <- comparison_br %>%
  group_by(Constituent) %>%
  summarize(n                       = n(),
            Mean_Log_Diff           = mean(Log_Diff),
            SD_Log_Diff             = sd(Log_Diff),
            Median_Pct_Diff         = median(Pct_Diff),
            Mean_network_distance_m = mean(network_distance_m),
            Max_network_distance_m  = max(network_distance_m),
            Mean_day_gap            = mean(day_gap),
            Max_day_gap             = max(day_gap),
            .groups = "drop")
print(comparison_summary_br)
print(comparison_br %>% filter(Constituent == "NPOC..mg.C.L.") %>% select(CampaignID, Site, nearest_predID, Conc_real, Conc_synth, Log_Diff), n=50)
print(comparison_br %>% filter(Constituent == "TDN..mg.N.L.") %>% select(CampaignID, Site, nearest_predID, Conc_real, Conc_synth, Log_Diff), n=50)
# Visual check: real vs. synthetic concentration, one point per site/campaign
ggplot(comparison_br, aes(x = Conc_real, y = Conc_synth, color = factor(CampaignID))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  #scale_x_log10() +
  #scale_y_log10() +
  facet_wrap(~ Constituent, scales = "free") +
  labs(x = "Real concentration", y = "Synthetic concentration",
       color = "Campaign", title = "br: real vs. synthetic concentration by site") +
  theme_bw()

## Results: I am satisfied that the synthetic data provides realistic values that should be directly comparable the real dataset. 


### Second, assess whether the pattern of variation holds true with synthetic data
# In the real data, CV nm > CV br for both DOC and TDN and both CVt and CVs.

comparison_nm$Catchment = "nm"
comparison_br$Catchment = "br"
synth_dat = bind_rows(comparison_br, comparison_nm)
synth_dat <- synth_dat %>%
  pivot_wider(
    names_from  = Constituent,
    values_from = c(Conc_real, Conc_synth)
  )
constituents_synth = c("Conc_real_NPOC..mg.C.L.", "Conc_synth_NPOC..mg.C.L.",
                       "Conc_real_TDN..mg.N.L.", "Conc_synth_TDN..mg.N.L.")

## calculate and compare CVt
CVt_synth_result <- temporal_cv(data = synth_dat,
                          concentration = constituents_synth,
                          site = "Site",
                          watershed = "Catchment",
                          digits = 2)
# One row per site × constituent
CVt_by_site.2 <- 
  CVt_synth_result$by_site %>%
  transmute(Catchment = Catchment,
            Constituent = constituent,
            Site,
            n_sites = n_used,
            CVt = temporal_cv)
print(CVt_by_site.2)
# summarize
CVt_by_site.3 <- CVt_by_site.2 %>%
  mutate(
    Source      = str_extract(Constituent, "real|synth"),
    Constituent = str_remove(Constituent, "^Conc_(real|synth)_")
  )
CVt_wide <- CVt_by_site.3 %>%
  select(Catchment, Constituent, Site, Source, CVt) %>%
  pivot_wider(names_from = Source, values_from = CVt, names_prefix = "CVt_")
CVt_summary <- CVt_wide %>%
  pivot_longer(cols = c(CVt_real, CVt_synth), names_to = "Source", values_to = "CVt") %>%
  mutate(Source = str_remove(Source, "^CVt_")) %>%
  group_by(Catchment, Constituent, Source) %>%
  summarise(
    n_sites  = sum(!is.na(CVt)),
    mean_CVt = mean(CVt, na.rm = TRUE),
    sd_CVt   = sd(CVt, na.rm = TRUE),
    .groups = "drop"
  )
CVt_summary
cvt_summary_npoc <- CVt_summary %>% filter(Constituent == "NPOC..mg.C.L.")
cvt_summary_tdn <- CVt_summary %>% filter(Constituent == "TDN..mg.N.L.")

## calculate and compare CVs
CVs_result <- spatial_cv(data = synth_dat,
                         concentration = constituents_synth,
                         event = "CampaignID",
                         digits = 2)
# One row per campaign × constituent
CVs_by_campaign.2 <- 
  CVs_result$by_event %>%
  transmute(#Catchment = Catchment,
            Constituent = constituent,
            CampaignID,
            n_sites = n_used,
            CVs = spatial_cv) %>%
  mutate(Catchment = str_split_i(CampaignID, pattern = "_", i = 1))
# summarize
CVs_by_campaign.3 <- CVs_by_campaign.2 %>%
  mutate(
    Source      = str_extract(Constituent, "real|synth"),
    Constituent = str_remove(Constituent, "^Conc_(real|synth)_")
  )
CVs_wide <- CVs_by_campaign.3 %>%
  select(Catchment, Constituent, CampaignID, Source, CVs) %>%
  pivot_wider(names_from = Source, values_from = CVs, names_prefix = "CVs_")
CVs_summary <- CVs_wide %>%
  pivot_longer(cols = c(CVs_real, CVs_synth), names_to = "Source", values_to = "CVs") %>%
  mutate(Source = str_remove(Source, "^CVs_")) %>%
  group_by(Catchment, Constituent, Source) %>%
  summarise(
    n_sites  = sum(!is.na(CVs)),
    mean_CVs = mean(CVs, na.rm = TRUE),
    sd_CVs   = sd(CVs, na.rm = TRUE),
    .groups = "drop"
  )
CVs_summary
CVs_summary_npoc <- CVs_summary %>% filter(Constituent == "NPOC..mg.C.L.")
CVs_summary_tdn <- CVs_summary %>% filter(Constituent == "TDN..mg.N.L.")

# plot
CVt_CVs_biplot_df <- CVt_summary %>%
  inner_join(CVs_summary, by = c("Catchment", "Constituent", "Source"))
ggplot(CVt_CVs_biplot_df, aes(x = mean_CVt, y = mean_CVs, color = Source, shape = Catchment)) +
  geom_errorbar(aes(ymin = mean_CVs - sd_CVs, ymax = mean_CVs + sd_CVs), width = 0) +
  geom_errorbarh(aes(xmin = mean_CVt - sd_CVt, xmax = mean_CVt + sd_CVt), height = 0) +
  geom_point(size = 3) +
  facet_wrap(~ Constituent) +
  labs(x = "Mean CVt", y = "Mean CVs") +
  theme_bw()

# result: synthetic data shows a similar pattern to real data (higher variance in nm, and higher CVs than CVt), except for TDN in br. For this one, the synthetic data has notably higher CVs than the real data and the sd of CVs is very small. 


### Third, combine real and synthetic data to total 20 campaigns.

extend_campaigns_nm <- 10:26
stopifnot(length(extend_campaigns_nm) == 17)
# --- Real portion: relabel each real CampaignNum onto its matched synthetic slot ---
real_portion_nm <- clean_nm %>%
  mutate(CampaignID_original = CampaignID) %>%
  left_join(
    campaign_match_nm %>%
      select(CampaignID_original = CampaignID, CampaignNum = nearest_CampaignNum),  # <-- adjust "RealCampaign" if named differently
    by = "CampaignID_original"
  ) %>%
  mutate(Source = "real") %>%
  select(Site, CampaignNum, CampaignID_original, Source, all_of(constituents))
# --- Synthetic portion: each Site's matched upstream predpt, campaigns 10-26 ---
extended_nm_wide <- extended_nm %>%
  select(predID, CampaignNum, Constituent, Conc) %>%
  pivot_wider(names_from = Constituent, values_from = Conc)
synth_portion_nm <- site_match_nm_upstream %>%
  select(Site, predID = nearest_predID) %>%
  inner_join(extended_nm_wide, by = "predID", relationship = "many-to-many") %>%
  filter(CampaignNum %in% extend_campaigns_nm) %>%
  mutate(CampaignID_original = NA_character_, Source = "synth") %>%
  select(Site, CampaignNum, CampaignID_original, Source, all_of(constituents))
sites_in_all_real_campaigns_nm <- clean_nm %>%
  distinct(Site, CampaignID) %>%
  count(Site) %>%
  filter(n == 3) %>%
  pull(Site)
synth_portion_nm <- synth_portion_nm %>%
  filter(Site %in% sites_in_all_real_campaigns_nm)
# --- Combine, add a convenience 1:20 sequential index ---
combined_nm <- bind_rows(real_portion_nm, synth_portion_nm) %>%
  arrange(CampaignNum, Site) %>%
  mutate(CampaignSeq = match(CampaignNum, sort(unique(CampaignNum))))
# --- checks ---
count(combined_nm, Source)
count(combined_nm, CampaignNum) %>% arrange(CampaignNum)
combined_nm %>% count(Site, CampaignNum) %>% filter(n > 1)
combined_nm %>% filter(if_any(all_of(constituents), is.na))

extend_campaigns_br <- 10:26
stopifnot(length(extend_campaigns_br) == 17)
# --- Real portion: relabel each real CampaignNum onto its matched synthetic slot ---
real_portion_br <- clean_br %>%
  mutate(CampaignID_original = CampaignID) %>%
  left_join(
    campaign_match_br %>%
      select(CampaignID_original = CampaignID, CampaignNum = nearest_CampaignNum),  # <-- adjust "RealCampaign" if named differently
    by = "CampaignID_original"
  ) %>%
  mutate(Source = "real") %>%
  select(Site, CampaignNum, CampaignID_original, Source, all_of(constituents))
# --- Synthetic portion: each Site's matched upstream predpt, campaigns 10-26 ---
extended_br_wide <- extended_br %>%
  select(predID, CampaignNum, Constituent, Conc) %>%
  pivot_wider(names_from = Constituent, values_from = Conc)
synth_portion_br <- site_match_br_upstream %>%
  select(Site, predID = nearest_predID) %>%
  inner_join(extended_br_wide, by = "predID", relationship = "many-to-many") %>%
  filter(CampaignNum %in% extend_campaigns_br) %>%
  mutate(CampaignID_original = NA_character_, Source = "synth") %>%
  select(Site, CampaignNum, CampaignID_original, Source, all_of(constituents))
sites_in_all_real_campaigns_br <- clean_br %>%
  distinct(Site, CampaignID) %>%
  count(Site) %>%
  filter(n == 3) %>%
  pull(Site)
synth_portion_br <- synth_portion_br %>%
  filter(Site %in% sites_in_all_real_campaigns_br)
# --- Combine, add a convenience 1:20 sequential index ---
combined_br <- bind_rows(real_portion_br, synth_portion_br) %>%
  arrange(CampaignNum, Site) %>%
  mutate(CampaignSeq = match(CampaignNum, sort(unique(CampaignNum))))
# --- checks ---
count(combined_br, Source)
count(combined_br, CampaignNum) %>% arrange(CampaignNum)
combined_br %>% count(Site, CampaignNum) %>% filter(n > 1)
combined_br %>% filter(if_any(all_of(constituents), is.na))

# combine catchments
combined_nm$Catchment = "nm"
combined_br$Catchment = "br"
combined = bind_rows(combined_br, combined_nm)

### Fourth, calculate CVt and CVs with ~20 campaigns.

CVt_combined_result <- temporal_cv(data = combined,
                                concentration = constituents,
                                site = "Site",
                                watershed = "Catchment",
                                digits = 2)
# One row per site × constituent
CVt_combined_result.2 <- 
  CVt_combined_result$by_site %>%
  transmute(Catchment = Catchment,
            Constituent = constituent,
            Site,
            n_sites = n_used,
            CVt = temporal_cv)
print(CVt_combined_result.2)

CVt_combined_result.3 = CVt_combined_result.2 %>% filter(n_sites>10)


CVs_result <- spatial_cv(data = combined,
                         concentration = constituents,
                         event = "CampaignNum",
                         watershed = "Catchment",
                         digits = 2)
# One row per campaign × constituent
CVs_by_campaign.2 <- 
  CVs_result$by_event %>%
  transmute(Catchment = Catchment,
    Constituent = constituent,
    CampaignNum,
    n_sites = n_used,
    CVs = spatial_cv) 


#### PART A.5 -- overall average for conceptual figure ####

# real data only
CVt_plot_data %>% group_by(Catchment, Constituent) %>% summarise(mean = mean(CVt), sd = sd(CVt))

# real + synthetic data
CVt_combined_result.3 %>% group_by(Catchment, Constituent) %>% summarise(mean = mean(CVt), sd = sd(CVt))
CVs_by_campaign.2 %>% group_by(Catchment, Constituent) %>% summarise(mean = mean(CVs), sd = sd(CVs))
