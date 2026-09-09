#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Lauren Giggy (with help adapting correlation code for individual constituents into a function to sweep through multiple constituents from ChatGPT version GPT-5.6 Luna 2026-07-27 and 2026-08-02 reviewed for consistency and edited by L. Giggy)
# (Also adapted from code authored by Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster))
# Last update: Lauren Giggy, 2026-09-09

# Requires: 02_build_synthetic_data.R must be run first (produces data/nm_clean.csv, data/nm_field_setups.rds, and data/nm_synthetic_extended.csv), and spatiotemporal_helpers.R must be in the same folder as this script.

# This script demonstrates how to calculate Spatial Stability, a metric of how consistent solute concentrations rankings are over a series of synoptic campaigns (one value per sampling event). 
# This script demonstrates how to calculate Spatial Stability with 2 different methods
# It also performs a sensitivity analysis to answer the question: how much does a Spatial Stability estimate vary depending on how many sites you sample? This analysis uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possible fits to a SSN2 spatial covariance model. 

# This script is organized as:
#   PART A -- SPpairs: computed directly from the real amd extended nm data
#   PART B -- SPfwmc: computed directly from the real amd extended nm data
#   PART C -- Monte Carlo sensitivity sweep
#   PART D -- Visualize results to compare methods

# Outputs:
#   data/SPpairs_by_campaign.csv:  real SPpairs, one row per real campaign x constituent
#   data/SPpairs_observed.csv: Mean_SPpairs and SD_SPpairs across those real campaigns
#   data/SPpairs_by_campaign_ext.csv:  extended SPpairs, one row per real campaign x constituent
#   data/SPpairs_observed_ext.csv: Mean_SPpairs and SD_SPpairs across extended campaigns
#   data/SPfwmc_by_campaign.csv:  real SPfwmc, one row per real campaign x constituent
#   data/SPfwmc_observed.csv: Mean_SPfwmc and SD_SPfwmc across those real campaigns
#   data/SPfwmc_by_campaign_ext.csv:  extended SPfwmc, one row per real campaign x constituent
#     - note Q data not available for SPfwmc_by_campaign_ext, so values are based on a mean concentration, not flow weighted
#   data/SPfwmc_observed_ext.csv: Mean_SPfwmc and SD_SPfwmc across extended campaigns


#   plots/SPpairs_SPfwmc_by_campaign.png: SPpairs method results and SPfwmc method results for observed data side by side
#   plots/SPpairs_SPfwmc_extended.png: SPpairs method results and SPfwmc method results for extended data side by side


#### Packages ####
library(dplyr)
library(googledrive)
library(lubridate)
library(Hmisc)
library(broom)
library(purrr)
library(tidyverse)
library(ggplot2)
library(cowplot)

source("spatiotemporal_helpers.R")













### Read in Data
# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM_BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toy = read.csv("drivedata/toy.csv")


### Read in Data
# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM_BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toy = read.csv("drivedata/toy.csv")

#Adjust datetime
#View(toy)
toy$Date <- lubridate::mdy(toy$Date)

#add area normalized streamflow values - Q must be in units of L/s and area in units of m2 for this formula
toy <- mutate(toy, q_mm_day = (Q* 86400)/Area.m2)

#Subset Toy Data into Projects 
NM_toy <- subset(toy, Project %in% c('nm'))
BR_toy <- subset(toy, Project %in% c('br'))

#Add event IDs to data based on date of data collection. 
#This allows samples collected within 2 days of each other to be combined into 1 event. 
#Change the "units" or "time dif" as necessary for study design. 
NM_toy <-  NM_toy %>%
  arrange(Date) %>%  # make sure data are ordered
  mutate(
    time_diff = as.numeric(difftime(Date, lag(Date), units = "days")),
    new_event = ifelse(is.na(time_diff) | time_diff > 2, 1, 0),
    event_id = cumsum(new_event)
  ) %>%
  dplyr::select(-time_diff, -new_event)

BR_toy <-  BR_toy %>%
  arrange(Date) %>%  # make sure data are ordered
  mutate(
    time_diff = as.numeric(difftime(Date, lag(Date), units = "days")),
    new_event = ifelse(is.na(time_diff) | time_diff > 2, 1, 0),
    event_id = cumsum(new_event)
  ) %>%
  dplyr::select(-time_diff, -new_event)

#Check for repeats
NM_toy |>
  summarise(n = dplyr::n(), .by = c(Site, event_id)) |>
  filter(n > 1)

BR_toy |>
  summarise(n = dplyr::n(), .by = c(Site, event_id)) |>
  filter(n > 1)


### median spearman rank and summary table function 
calc_event_cor <- function(data, analyte, min_samples = 8) { 
  #adjust min_samples as preferred
  
  # Create Solute × Event matrix and set up for correlation calculation
  mat <- data %>%
    dplyr::select(Site, event_id, all_of(analyte)) %>% #keep only the columns needed for calculation
    drop_na() %>% #remove NAs
    pivot_wider(names_from = event_id, values_from = all_of(analyte)) %>% #rotate the dataframe so each column is a synoptic event
    arrange(Site) %>% #arrange matching Sites
    select(-Site) %>% #removes Site for correlation 
    as.matrix()
  
  # Calculate the number of paired observations for each event
  n_mat <- outer(
    seq_len(ncol(mat)),
    seq_len(ncol(mat)),
    Vectorize(function(i, j) {
      sum(complete.cases(mat[, c(i, j)]))
    }))
  
  # Calculate event-by-event Spearman correlations
  cor_mat <- cor(mat, method = "spearman",
                 use = "pairwise.complete.obs")
  
  # Remove correlations based on too few paired observations
  cor_mat[n_mat < min_samples] <- NA
  
  #Remove event correlations with itself 
  diag(cor_mat) <- NA
  
  #create a table of the median spearman rank, MARGIN = 2 means it will apply across columns (= 1 applies across rows)
  tibble(event_id = as.numeric(colnames(cor_mat)),
         median_spearman = apply(cor_mat, MARGIN = 2, median, na.rm = TRUE), analyte = analyte)} 

#define analytes
analytes <- c(
  "NPOC..mg.C.L.",
  "TDN..mg.N.L.",
  "q_mm_day")

#Run function on data
NM_toy_sp_pairs <- map_dfr(analytes, ~calc_event_cor(NM_toy, .x))
BR_toy_sp_pairs <- map_dfr(analytes, ~calc_event_cor(BR_toy, .x))


#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

n_iter    <- 10000                                  # Monte Carlo iterations
site_grid <- c(3:50, 75, 100, 150, 300)              # site counts to test

actual_n_sites <- 23  # real number of nm sites sampled -- update if that changes

set.seed(42)

clean            <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
field_setups     <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)


#### PART A -- Calculate SpSt pairs of real toy and extended dataset ####
# Each site's concentration averaged across its available real campaigns, then CV taken across sites.

Constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.") #name constituents

calc_SPpairs <- function(data, Constituent, min_samples = 0) {
  
  # Create Site × Event matrix
  mat <- data %>%
    dplyr::select(Site, CampaignID, all_of(Constituent)) %>%
    drop_na() %>%
    pivot_wider(
      names_from = CampaignID,
      values_from = all_of(Constituent)
    ) %>%
    arrange(Site) %>%
    dplyr::select(-Site) %>%
    as.matrix()
  
  
  # Number of paired observations for each event pair
  n_mat <- outer(
    seq_len(ncol(mat)),
    seq_len(ncol(mat)),
    Vectorize(function(i, j) {
      sum(complete.cases(mat[, c(i, j)]))
    })
  )
  
  # Event-by-event Spearman correlations
  cor_mat <- cor(
    mat,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  # Remove correlations based on too few paired observations
  cor_mat[n_mat < min_samples] <- NA
  
  diag(cor_mat) <- NA
  
  #create a table of the median spearman rank, MARGIN = 2 means it will apply across columns (= 1 applies across rows)
  tibble(
    CampaignID = (colnames(cor_mat)),
    median_spearman = apply(cor_mat, MARGIN = 2, median, na.rm = TRUE),
    Constituent = Constituent
  )
} #this is a function that calculate the SPpairs for each consituent

SPpairs_by_campaign <- map_dfr(Constituents, ~calc_SPpairs(clean, .x)) #apply function to toy data

print(SPpairs_by_campaign) #check results

SPpairs_by_campaign <- SPpairs_by_campaign  %>% rename(SPpairs = median_spearman) #rename results column for clarity

write_csv(SPpairs_by_campaign, file.path(data_out_dir, "SPpairs_by_campaign.csv"))

## Summarized across campaigns per constituent
SPpairs_observed <- SPpairs_by_campaign %>%
  filter(!is.na(SPpairs)) %>%
  group_by(Constituent) %>%
  dplyr::summarize(n_campaigns = n(), Mean_SPpairs = mean(SPpairs),
                   SD_SPpairs = if (n() >= 2) sd(SPpairs) else NA_real_, .groups = "drop")

print(SPpairs_observed)#check results

write_csv(SPpairs_observed, file.path(data_out_dir, "SPpairs_observed.csv"))

###  calcululate SPpairs on extended data   ###

#pivot df to match structure of toy dataset
synthetic_extended_wide <-synthetic_extended %>% tidyr::pivot_wider( 
  names_from = Constituent,
  values_from = Conc) 

Constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.") #name constituents

calc_SPpairs_ext <- function(data, Constituent, min_samples = 0) {
  
  # Create predID × Event matrix
  mat <- data %>%
    dplyr::select(predID, CampaignNum, all_of(Constituent)) %>%
    drop_na() %>%
    pivot_wider(
      names_from = CampaignNum,
      values_from = all_of(Constituent)
    ) %>%
    arrange(predID) %>%
    dplyr::select(-predID) %>%
    as.matrix()
  
  
  # Number of paired observations for each event pair
  n_mat <- outer(
    seq_len(ncol(mat)),
    seq_len(ncol(mat)),
    Vectorize(function(i, j) {
      sum(complete.cases(mat[, c(i, j)]))
    })
  )
  
  # Event-by-event Spearman correlations
  cor_mat <- cor(
    mat,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  # Remove correlations based on too few paired observations
  cor_mat[n_mat < min_samples] <- NA
  
  diag(cor_mat) <- NA
  
  #create a table of the median spearman rank, MARGIN = 2 means it will apply across columns (= 1 applies across rows)
  tibble(
    CampaignNum = (colnames(cor_mat)),
    median_spearman = apply(cor_mat, MARGIN = 2, median, na.rm = TRUE),
    Constituent = Constituent
  )
}  # SPpairs function adjusted to match column names of extended data 

SPpairs_by_campaign_ext <- map_dfr(Constituents, ~calc_SPpairs_ext(synthetic_extended_wide, .x)) #apply function to extended data

print(SPpairs_by_campaign_ext) #check data

SPpairs_by_campaign_ext <- SPpairs_by_campaign_ext  %>% rename(SPpairs = median_spearman) #rename results column for clarity

write_csv(SPpairs_by_campaign_ext, file.path(data_out_dir, "SPpairs_by_campaign_ext.csv"))

## Summarized across campaigns per constituent
SPpairs_observed_ext <- SPpairs_by_campaign_ext %>%
  filter(!is.na(SPpairs)) %>%
  group_by(Constituent) %>%
  dplyr::summarize(n_campaigns = n(), Mean_SPpairs = mean(SPpairs),
                   SD_SPpairs = if (n() >= 2) sd(SPpairs) else NA_real_, .groups = "drop")

print(SPpairs_observed_ext)#check results

write_csv(SPpairs_observed_ext, file.path(data_out_dir, "SPpairs_observed_ext.csv"))


#### PART B -- Calculate SPfwmc of real toy and extended dataset ####
# Each site's concentration averaged across its available real campaigns, then CV taken across sites.

Constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.") #name constituents

calc_SPfwmc_event_cor <- function(data, Constituent) {
  
  # Calculate flow-weighted mean by Site
  data_fwmc <- data %>%
    group_by(Site) %>%
    summarise(
      FWMC = weighted.mean(
        .data[[Constituent]],
        w = Q,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # Put CampaignID values into columns
  data_wide <- data %>%
    dplyr::select(Site, CampaignID, value = all_of(Constituent)) %>%
    drop_na() %>%
    pivot_wider(
      names_from = CampaignID,
      values_from = value
    ) %>%
    arrange(Site)
  
  # Merge FWMC with campaign data
  data_fwmc_merge <- merge(
    data_fwmc,
    data_wide,
    by = "Site"
  )
  
  data_fwmc_merge$Site <- NULL
  
  data_fwmc_merge[] <- lapply(
    data_fwmc_merge,
    as.numeric
  )
  
  # Spearman correlation
  cor_mat <- cor(
    data_fwmc_merge,
    method = "spearman",
    use = "pairwise.complete.obs")
  
  # Don't correlate FWMC with itself
  diag(cor_mat) <- NA
  
  # Extract correlations with FWMC
  cor_df <- as.data.frame(cor_mat)
  cor_df_final <- cor_df[, "FWMC", drop = FALSE]
  cor_df_final$CampaignID <- rownames(cor_df_final)
  cor_df_final <- na.omit(cor_df_final)
  
  # Add constituent name
  cor_df_final$Constituent <- Constituent
  
  cor_df_final
} #this is a function that calculate the SPfwmc for each consituent\

SPfwmc_by_campaign <- map_dfr(Constituents, ~calc_SPfwmc_event_cor(clean, .x)) #apply function to toy data

print(SPfwmc_by_campaign) #check results

SPfwmc_by_campaign <- SPfwmc_by_campaign  %>% rename(SPfwmc = FWMC) #rename results column for clarity

write_csv(SPfwmc_by_campaign, file.path(data_out_dir, "SPfwmc_by_campaign.csv"))

###  calcululate SPpairs on extended data   ###

#pivot df to match structure of toy dataset
synthetic_extended_wide <-synthetic_extended %>% tidyr::pivot_wider(
  names_from = Constituent, values_from = Conc)

Constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

# note: Q data for extended data not available, taking mean C not flow weighted
calc_SPfwmc_event_cor_ext <- function(data, Constituent) {
  
  # Calculate flow-weighted mean by predID
  data_fwmc <- data %>%
    group_by(predID) %>%
    summarise(
      FWMC = mean(
        .data[[Constituent]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # Put CampaignNum values into columns
  data_wide <- data %>%
    dplyr::select(predID, CampaignNum, value = all_of(Constituent)) %>%
    drop_na() %>%
    pivot_wider(
      names_from = CampaignNum,
      values_from = value
    ) %>%
    arrange(predID)
  
  # Merge FWMC with campaign data
  data_fwmc_merge <- merge(
    data_fwmc,
    data_wide,
    by = "predID"
  )
  
  data_fwmc_merge$predID <- NULL
  
  data_fwmc_merge[] <- lapply(
    data_fwmc_merge,
    as.numeric
  )
  
  # Spearman correlation
  cor_mat <- cor(
    data_fwmc_merge,
    method = "spearman",
    use = "pairwise.complete.obs")
  
  # Don't correlate FWMC with itself
  diag(cor_mat) <- NA
  
  # Extract correlations with FWMC
  cor_df <- as.data.frame(cor_mat)
  cor_df_final <- cor_df[, "FWMC", drop = FALSE]
  cor_df_final$CampaignNum <- rownames(cor_df_final)
  cor_df_final <- na.omit(cor_df_final)
  
  # Add constituent name
  cor_df_final$Constituent <- Constituent
  
  cor_df_final
}  # SPfwmc function adjusted to match column names of extended data 

SPfwmc_by_campaign_ext <- map_dfr(Constituents, ~calc_SPfwmc_event_cor_ext(synthetic_extended_wide, .x))  #apply function to extended data

print(SPfwmc_by_campaign_ext) #check results

SPfwmc_by_campaign_ext <- SPfwmc_by_campaign_ext  %>% rename(SPfwmc = FWMC) #rename results column for clarity

write_csv(SPfwmc_by_campaign_ext, file.path(data_out_dir, "SPfwmc_by_campaign.csv"))



#### PART C -- Monte Carlo spatial sensitivity analysis ####

spatial_mc_results <- list()

for (cc in names(field_setups)) {
  fs <- field_setups[[cc]]
  m <- fs$setup$m
  this_site_grid <- site_grid[site_grid <= m]
  if (length(this_site_grid) == 0) {
    cat("  SKIP", cc, "-- every site_grid value exceeds m =", m, "predpts.\n")
    next
  }
  
  cat("\n---", cc, "[spatial MC,", n_iter, "iterations, up to", max(this_site_grid), "sites] ---\n")
  
  ## Calc reference CVs
  ref_cvs <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(CampaignNum) %>%
    summarize(CVs = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVs) %>% mean()
  
  ## One fresh, spatially-correlated noise draw per Monte Carlo iteration; drawn all at once (an m x n_iter matrix in one call). t_days is just a placeholder vector of the right length.
  noise_pool <- draw_field(fs$setup, t_days = seq_len(n_iter), rho = 0)
  sim_mat <- exp(fs$setup$mu_cond + noise_pool)  # m x n_iter, natural scale
  
  for (n in this_site_grid) {
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx <- sample(m, n)
      vals <- sim_mat[idx, i]
      cvs_draws[i] <- sd(vals) / mean(vals)
    }
    spatial_mc_results[[paste(cc, n)]] <- tibble(Constituent = cc, N = n, CVs = cvs_draws, Ref_CVs = ref_cvs)
  }
}

spatial_mc_results <- bind_rows(spatial_mc_results)

# present and save results
spatial_mc_summary <- spatial_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs = first(Ref_CVs), .groups = "drop") %>%
  mutate(Pct_Bias = (Median_CVs - Ref_CVs) / Ref_CVs * 100)


print(spatial_mc_summary, n = Inf)

#write_csv(spatial_mc_summary, file.path(data_out_dir, "CVs_MC_summary.csv"))

## Creates one row per constituent
bias_annot <- spatial_mc_summary %>%
  filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.1f%% bias", N, Pct_Bias))

## First plot: how CVs changes with sample size, with % bias of actual sample size annotated
p_cvs <- ggplot(spatial_mc_summary, aes(x = N, y = Median_CVs)) +
  geom_ribbon(aes(ymin = P05_CVs, ymax = P95_CVs), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVs), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dashed") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites)", y = "CVs (median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity",
       subtitle = "Blue = reference CVs; red = actual sample size") +
  theme_bw()
#ggsave(file.path(plot_dir, "CVs_MC_plot.png"), p_cvs, width = 8, height = 6, dpi = 150)

## Second plot: percent bias and precision (SD across Monte Carlo iterations), both as a function of N. Pct_Bias shows how far off the median estimate tends to be; SD_CVs shows how much that estimate itself changes at a given sample size. Pivoted to long format so both metrics can share one facet grid (metric x constituent) rather than needing a dual y-axis.
sensitivity_long <- spatial_mc_summary %>%
  select(Constituent, N, Pct_Bias, SD_CVs) %>%
  pivot_longer(c(Pct_Bias, SD_CVs), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric, Pct_Bias = "Percent bias vs Ref_CVs (%)",
                         SD_CVs = "SD of CVs across MC iterations"))
p_sensitivity <- ggplot(sensitivity_long, aes(x = N, y = Value)) +
  geom_line() + geom_point(size = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dashed") +
  facet_grid(Metric ~ Constituent, scales = "free_y") +
  labs(x = "Sample size (N sites)", y = NULL,
       title = "Spatial sensitivity: bias and precision vs. sample size",
       subtitle = "Red = actual sample size") +
  theme_bw()





#### PART D - Visualize ------------------------------------------------------

###   Compare SPpairs and SPfwmc results for Toy data  ###
SPpairs<- ggplot(data = SPpairs_by_campaign, aes(x = CampaignID, y = SPpairs, color = Constituent))+
  geom_point(size =2)+
  ylim(-1,1)+
  labs(title = "SPpairs Toy Data")+
  theme(legend.position="right")+
  theme_classic(base_size =10) 

SPfwmc <- ggplot(data = SPfwmc_by_campaign, aes(x = CampaignID, y = SPfwmc, color = Constituent))+
  geom_point(size =2)+
  ylim(-1,1)+
  labs(title = "SPfwmc Toy Data")+
  theme(legend.position="right") +
  theme_classic(base_size =10 )

SPpairs_SPfwmc_by_campaign <- plot_grid(SPpairs, SPfwmc)

plot(SPpairs_SPfwmc_by_campaign) #View side by side plots

ggsave(file.path(plot_dir, "SPpairs_SPfwmc_by_campaign.png"), SPpairs_SPfwmc_by_campaign, width = 8, height = 4, dpi = 150)


###   Compare SPpairs and SPfwmc results for extended data  ###
SPpairs_ext <- ggplot(data = SPpairs_by_campaign_ext, aes(x = CampaignNum, y = SPpairs, color = Constituent))+
  geom_point(size =2)+
  ylim(-1,1)+
  labs(title = "SPpairs Extended Data")+
  theme(legend.position="right") +
  theme_classic(base_size =10 )

SPfwmc_ext <- ggplot(data = SPfwmc_by_campaign_ext, aes(x = CampaignNum, y = SPfwmc, color = Constituent))+
  geom_point(size =2)+
  ylim(-1,1)+
  labs(title = "SPfwmc Extended Data")+
  theme(legend.position="right") +
  theme_classic(base_size =10 )

SPpairs_SPfwmc_by_campaign <- plot_grid(SPpairs_ext, SPfwmc_ext)

plot(SPpairs_SPfwmc_by_campaign) #View side by side plots

ggsave(file.path(plot_dir, "SPpairs_SPfwmc_extended.png"), SPpairs_SPfwmc_by_campaign, width = 8, height = 4, dpi = 150)



