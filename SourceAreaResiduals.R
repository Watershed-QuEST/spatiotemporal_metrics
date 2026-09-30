#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): William Mejia, 29-Sep-2026

# This script calculates source area residual (SAR) for QuEST toy
# NM-BR dataset, used for the planned commentary manuscript on
# spatiotemporal metrics for assessing variance in stream chemistry
# across stream networks.

##### Load packages #####

library(tidyverse)
library(googledrive)
library(lubridate)
library(dataRetrieval)
library(ggpmisc)


##### Read in data #####

drive_folder <- "https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"
toy_files    <- drive_ls(drive_get(drive_folder))

drive_download(
  file      = toy_files$id[toy_files$name == "NM-BR Toy dataset.csv"],
  path      = "drivedata/toy.csv",
  overwrite = TRUE
)

toy <- read.csv("drivedata/toy.csv")


#===============================================================================
# DATA PREPARATION
#===============================================================================

# Instantaneous mass loads
# Mass = Q * concentration

toy <- toy %>%
  mutate(
    Project = str_to_upper(Project),    # can we have them be uppercase in the og data? lol
    npoc_mass = Q * NPOC..mg.C.L.,
    tdn_mass  = Q * TDN..mg.N.L.
  )


##### Solute columns #####

mass_cols <- c("npoc_mass", "tdn_mass")


#===============================================================================
# SITE-LEVEL SUMMARY
#===============================================================================
#
# The regression models are fit using site-level mean values.
#
# For each site:
#   - Project = watershed
#   - Area = site contributing area (upstream of sampling point)
#   - Q = mean discharge
#   - mass = mean mass load
#
# The AQ model therefore uses mean site discharge when fitting the regression.
#===============================================================================

site_summary <- toy %>%
  group_by(Site) %>%
  summarise(
    
    Project = first(Project),
    area_m2 = first(Area.m2),
    
    across(
      all_of(c("Q", mass_cols)),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    
    # Site-level log-transformed predictors
    logArea = log10(area_m2),
    logQ    = log10(ifelse(Q_mean > 0, Q_mean, NA)),
    
    # Site-level log-transformed mass
    logNPOC = log10(ifelse(npoc_mass_mean > 0, npoc_mass_mean, NA)),
    logTDN  = log10(ifelse(tdn_mass_mean > 0, tdn_mass_mean, NA))
  )


#===============================================================================
# FULL DATASET LOG TRANSFORMS
#===============================================================================
#
# These are used for calculating observation-level residuals.
#
# Area is constant within each site, while Q can vary among observations.
#===============================================================================

toy_log <- toy %>%
  mutate(
    
    logArea = log10(Area.m2),
    
    logQ = log10(
      ifelse(Q > 0, Q, NA)
    ),
    
    logNPOC = log10(
      ifelse(npoc_mass > 0, npoc_mass, NA)
    ),
    
    logTDN = log10(
      ifelse(tdn_mass > 0, tdn_mass, NA)
    )
  )


#===============================================================================
# AREA-ONLY (A) MODELS
#===============================================================================
#
# Model formulation:
#
#   log10(M_i) = alpha_A + beta_A * log10(A_i) + epsilon_i
#
# The fitted prediction is:
#
#   yhat_i,A = alphahat_A + betahat_A * log10(A_i)
#
# SAR:
#
#   SAR_ij,A = log10(M_ij) - yhat_i,A
#
#===============================================================================

lm_A <- list(
  
  logNPOC = lm(
    logNPOC ~ logArea,
    data = site_summary
  ),
  
  logTDN = lm(
    logTDN ~ logArea,
    data = site_summary
  )
)


#===============================================================================
# AREA + DISCHARGE (AQ) MODELS (this is the model the NM team uses IRL)
#===============================================================================
#
# Model formulation:
#
#   log10(M_i) =
#       alpha_AQ
#       + beta_A * log10(A_i)
#       + beta_Q * log10(Q_i)
#       + epsilon_i
#
# The fitted prediction is:
#
#   yhat_i,AQ =
#       alphahat_AQ
#       + betahat_A * log10(A_i)
#       + betahat_Q * log10(Q_i)
#
# For observation-level predictions, individual observation Q is used.
#
# SAR:
#
#   SAR_ij,AQ = log10(M_ij) - yhat_ij,AQ
#
#===============================================================================

lm_AQ <- list(
  
  logNPOC = lm(
    logNPOC ~ logArea + logQ,
    data = site_summary
  ),
  
  logTDN = lm(
    logTDN ~ logArea + logQ,
    data = site_summary
  )
)


#===============================================================================
# PREDICTIONS + RESIDUALS
#===============================================================================

#-------------------------------------------------------------------------------
# Area-only predictions and SARs
#-------------------------------------------------------------------------------

for (solute in names(lm_A)) {
  
  pred_name <- paste0("PL_", solute, "_A")
  res_name  <- paste0("SAR_", solute, "_A")
  
  toy_log[[pred_name]] <-
    predict(
      lm_A[[solute]],
      newdata = toy_log
    )
  
  toy_log[[res_name]] <-
    toy_log[[solute]] - toy_log[[pred_name]]
}


#-------------------------------------------------------------------------------
# Area + discharge predictions and SARs
#-------------------------------------------------------------------------------

for (solute in names(lm_AQ)) {
  
  pred_name <- paste0("PL_", solute, "_AQ")
  res_name  <- paste0("SAR_", solute, "_AQ")
  
  toy_log[[pred_name]] <-
    predict(
      lm_AQ[[solute]],
      newdata = toy_log
    )
  
  toy_log[[res_name]] <-
    toy_log[[solute]] - toy_log[[pred_name]]
}


#===============================================================================
# View OBSERVATION-LEVEL RESIDUALS
#===============================================================================

obs_level_data <- toy_log %>%
  select(
    Site,
    Project,
    Area.m2,
    Q,
    logArea,
    logQ,
    
    # SARs
    starts_with("SAR_"),
    
    # Predictions
    starts_with("PL_"))

View(obs_level_data) #might be easier to separate model output

#===============================================================================
# RESIDUAL SUMMARIZATION
#===============================================================================


sar_cols <- grep(
  "^SAR_",
  names(toy_log),
  value = TRUE
)

residual_summary <- toy_log %>%
  group_by(
    Site,
    Project,
    Area.m2
  ) %>%
  summarise(
    
    across(
      all_of(sar_cols),
      
      list(
        
        mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        
        se = ~ sd(
          .x,
          na.rm = TRUE
        ) /
          sqrt(
            sum(!is.na(.x))
          )
      ),
      
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  )

View(residual_summary)

#===============================================================================
# PLOTS: AREA-SCALING RELATIONSHIPS BY WATERSHED
#===============================================================================

watersheds <- c("NM", "BR")

for (w in watersheds) {
  
  p <- ggplot(
    site_summary[
      site_summary$Project == w,
    ],
    aes(
      x = logArea,
      y = logNPOC
    )
  ) +
    
    geom_point(
      size = 3
    ) +
    
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black"
    ) +
    
    stat_poly_eq(
      aes(
        label = paste(
          ..eq.label..,
          ..rr.label..,
          sep = "~~~"
        )
      ),
      formula = y ~ x,
      parse = TRUE,
      size = 5
    ) +
    
    labs(
      x = expression(
        Log[10] ~ Subcatchment ~ Area ~ (m^2)
      ),
      
      y = expression(
        Log[10] ~ Site ~ Avg. ~ DOC ~ Mass
      ),
      
      title = paste0(
        w,
        ": Log-Log Relationship Between Site-Averaged DOC Mass\n",
        "and Subcatchment Area"
      )
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  print(p)
  
  ggsave(
    filename = paste0(
      "04_figures/",
      w,
      "_NPOC_area_scaling.png"
    ),
    
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
}

#must do the same for TDN

#===============================================================================
# PLOTS: RESIDUAL BAR PLOTS BY WATERSHED (this graph is confusing to me but og script has it, suggest change?)
#===============================================================================

ylims <- list(
  NM = c(-3, 1),
  BR = c(-3, 1.2)
)


for (w in watersheds) {
  
  #------------------------------------------------------------------------
  # Area-only SAR
  #------------------------------------------------------------------------
  
  p_A <- ggplot(
    residual_summary[
      residual_summary$Project == w,
    ]
  ) +
    
    geom_bar(
      aes(
        x = reorder(
          Site,
          SAR_logTDN_A_mean
        ),
        y = SAR_logTDN_A_mean
      ),
      
      stat = "identity",
      fill = "skyblue",
      alpha = 0.7
    ) +
    
    geom_errorbar(
      aes(
        x = Site,
        ymin = SAR_logTDN_A_mean -
          SAR_logTDN_A_se,
        ymax = SAR_logTDN_A_mean +
          SAR_logTDN_A_se
      ),
      
      width = 0.4,
      colour = "orange",
      alpha = 0.9,
      linewidth = 1.3
    ) +
    
    ylim(
      ylims[[w]]
    ) +
    
    labs(
      x = "Site",
      y = "TDN SAR (Area only)",
      title = paste0(
        w,
        ": Site-Level TDN SARs - Area Only"
      )
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  print(p_A)
  
  ggsave(
    filename = paste0(
      "04_figures/",
      w,
      "_TDN_SAR_A.png"
    ),
    
    plot = p_A,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  
  #------------------------------------------------------------------------
  # Area + discharge SAR
  #------------------------------------------------------------------------
  
  p_AQ <- ggplot(
    residual_summary[
      residual_summary$Project == w,
    ]
  ) +
    
    geom_bar(
      aes(
        x = reorder(
          Site,
          SAR_logTDN_AQ_mean
        ),
        y = SAR_logTDN_AQ_mean
      ),
      
      stat = "identity",
      fill = "skyblue",
      alpha = 0.7
    ) +
    
    geom_errorbar(
      aes(
        x = Site,
        ymin = SAR_logTDN_AQ_mean -
          SAR_logTDN_AQ_se,
        ymax = SAR_logTDN_AQ_mean +
          SAR_logTDN_AQ_se
      ),
      
      width = 0.4,
      colour = "orange",
      alpha = 0.9,
      linewidth = 1.3
    ) +
    
    ylim(
      ylims[[w]]
    ) +
    
    labs(
      x = "Site",
      y = "TDN SAR (Area + discharge)",
      title = paste0(
        w,
        ": Site-Level TDN SARs - Area + Discharge"
      )
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  print(p_AQ)
  
  ggsave(
    filename = paste0(
      "04_figures/",
      w,
      "_TDN_SAR_AQ.png"
    ),
    
    plot = p_AQ,
    width = 6,
    height = 4,
    dpi = 300
  )
}

# must look into what a negative residuals represents considering we 
# have yet to see it in our real data. is it as common as portayed here?
