#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): William Mejia, 01-OCT-2026

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

# Calculate the mean of log-transformed mass loads and discharge
# for each site. These site-level means are used to fit the regressions.
# Note: I was previously averaging site values, then log-transforming.


site_summary <- toy %>%
  mutate(
    logQ = log10(ifelse(Q > 0, Q, NA_real_)),
    logNPOC = log10(ifelse(npoc_mass > 0, npoc_mass, NA_real_)),
    logTDN = log10(ifelse(tdn_mass > 0, tdn_mass, NA_real_))
  ) %>%
  group_by(Site) %>%
  summarise(
    Project = first(Project),
    area_m2 = first(Area.m2),

    logQ = if (all(is.na(logQ))) {
      NA_real_
    } else {
      mean(logQ, na.rm = TRUE)
    },

    logNPOC = if (all(is.na(logNPOC))) {
      NA_real_
    } else {
      mean(logNPOC, na.rm = TRUE)
    },

    logTDN = if (all(is.na(logTDN))) {
      NA_real_
    } else {
      mean(logTDN, na.rm = TRUE)
    },

    .groups = "drop"
  ) %>%
  mutate(
    logArea = log10(area_m2)
  )


#===============================================================================
# FULL DATASET LOG TRANSFORMS
#===============================================================================
#
# Individual observations must be log-transformed to calculate residuals.
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
# FIT SEPARATE MODELS FOR NM AND BR
#===============================================================================

# Previous code had model that would use data from both catchments, this creates
# individual models for each catchment. This in turn now creates 8 separate   
# models, 4 for each catchment, 2 for each solute (NPOC and TDN), 
# and 2 for each model type (Area-only and Area + Discharge).

watersheds <- c("NM", "BR")

# Area-only models
lm_A <- list()

# Area + discharge models
lm_AQ <- list()

for (w in watersheds) {
  
  project_data <- site_summary %>%
    filter(Project == w)
  
  # Area-only models
  lm_A[[w]] <- list(
    
    logNPOC = lm(
      logNPOC ~ logArea,
      data = project_data
    ),
    
    logTDN = lm(
      logTDN ~ logArea,
      data = project_data
    )
  )
  
  # Area + discharge models
  lm_AQ[[w]] <- list(
    
    logNPOC = lm(
      logNPOC ~ logArea + logQ,
      data = project_data
    ),
    
    logTDN = lm(
      logTDN ~ logArea + logQ,
      data = project_data
    )
  )
}


# Inspect model summaries
lm_A$NM$logNPOC |> summary()
lm_A$NM$logTDN  |> summary()
lm_A$BR$logNPOC |> summary()
lm_A$BR$logTDN  |> summary()

lm_AQ$NM$logNPOC |> summary()
lm_AQ$NM$logTDN  |> summary()
lm_AQ$BR$logNPOC |> summary()
lm_AQ$BR$logTDN  |> summary()



#===============================================================================
# PREDICTIONS + RESIDUALS BY PROJECT
#===============================================================================

for (w in watersheds) {
  
  project_rows <- which(toy_log$Project == w)
  
  project_observations <- toy_log[project_rows, ]
  
  #---------------------------------------------------------------------------
  # Area-only predictions and SARs
  #---------------------------------------------------------------------------
  
  for (solute in names(lm_A[[w]])) {
    
    pred_name <- paste0("PL_", solute, "_A_", w)
    res_name  <- paste0("SAR_", solute, "_A_", w)
    
    predictions <- predict(
      lm_A[[w]][[solute]],
      newdata = project_observations
    )
    
    toy_log[[pred_name]] <- NA_real_
    toy_log[[pred_name]][project_rows] <- predictions
    
    toy_log[[res_name]] <- NA_real_
    toy_log[[res_name]][project_rows] <-
      project_observations[[solute]] - predictions
  }
  
  
  #---------------------------------------------------------------------------
  # Area + discharge predictions and SARs
  #---------------------------------------------------------------------------
  
  for (solute in names(lm_AQ[[w]])) {
    
    pred_name <- paste0("PL_", solute, "_AQ_", w)
    res_name  <- paste0("SAR_", solute, "_AQ_", w)
    
    predictions <- predict(
      lm_AQ[[w]][[solute]],
      newdata = project_observations
    )
    
    toy_log[[pred_name]] <- NA_real_
    toy_log[[pred_name]][project_rows] <- predictions
    
    toy_log[[res_name]] <- NA_real_
    toy_log[[res_name]][project_rows] <-
      project_observations[[solute]] - predictions
  }
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
    starts_with("SAR_"),
    starts_with("PL_")
  )

# NM observations and NM model columns
View(
  obs_level_data %>%
    filter(Project == "NM") %>%
    select(
      Site, Project, Area.m2, Q, logArea, logQ,
      ends_with("_NM")
    )
)

# BR observations and BR model columns
View(
  obs_level_data %>%
    filter(Project == "BR") %>%
    select(
      Site, Project, Area.m2, Q, logArea, logQ,
      ends_with("_BR")
    )
)

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
# PLOTS: RESIDUAL BAR PLOTS BY WATERSHED
#===============================================================================

ylims <- list(
  NM = c(-3, 1),
  BR = c(-3, 1.2)
)

for (w in watersheds) {
  
  # Column names for this project's area-only SARs
  sar_A_mean <- paste0("SAR_logTDN_A_", w, "_mean")
  sar_A_se   <- paste0("SAR_logTDN_A_", w, "_se")
  
  # Column names for this project's area + discharge SARs
  sar_AQ_mean <- paste0("SAR_logTDN_AQ_", w, "_mean")
  sar_AQ_se   <- paste0("SAR_logTDN_AQ_", w, "_se")
  
  
  #---------------------------------------------------------------------------
  # Area-only SAR
  #---------------------------------------------------------------------------
  
  p_A <- ggplot(
    residual_summary %>%
      filter(Project == w)
  ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black"
    ) +
    
    geom_col(
      aes(
        x = reorder(Site, .data[[sar_A_mean]]),
        y = .data[[sar_A_mean]]
      ),
      fill = "skyblue",
      alpha = 0.7
    ) +
    
    geom_errorbar(
      aes(
        x = reorder(Site, .data[[sar_A_mean]]),
        ymin = .data[[sar_A_mean]] - .data[[sar_A_se]],
        ymax = .data[[sar_A_mean]] + .data[[sar_A_se]]
      ),
      width = 0.4,
      colour = "orange",
      alpha = 0.9,
      linewidth = 1.1
    ) +
    
    coord_cartesian(
      ylim = ylims[[w]]
    ) +
    
    labs(
      x = "Site",
      y = expression("Mean TDN SAR (" * log[10] * " scale)"),
      title = paste0(w, ": Site-Level TDN SARs — Area Only")
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  print(p_A)
  
  ggsave(
    filename = paste0(
      "04_figures/",
      w,
      "_TDN_SAR_A.png"
    ),
    plot = p_A,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  
  #---------------------------------------------------------------------------
  # Area + discharge SAR
  #---------------------------------------------------------------------------
  
  p_AQ <- ggplot(
    residual_summary %>%
      filter(Project == w)
  ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black"
    ) +
    
    geom_col(
      aes(
        x = reorder(Site, .data[[sar_AQ_mean]]),
        y = .data[[sar_AQ_mean]]
      ),
      fill = "skyblue",
      alpha = 0.7
    ) +
    
    geom_errorbar(
      aes(
        x = reorder(Site, .data[[sar_AQ_mean]]),
        ymin = .data[[sar_AQ_mean]] - .data[[sar_AQ_se]],
        ymax = .data[[sar_AQ_mean]] + .data[[sar_AQ_se]]
      ),
      width = 0.4,
      colour = "orange",
      alpha = 0.9,
      linewidth = 1.1
    ) +
    
    coord_cartesian(
      ylim = ylims[[w]]
    ) +
    
    labs(
      x = "Site",
      y = expression("Mean TDN SAR (" * log[10] * " scale)"),
      title = paste0(w, ": Site-Level TDN SARs — Area + Discharge")
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  print(p_AQ)
  
  ggsave(
    filename = paste0(
      "04_figures/",
      w,
      "_TDN_SAR_AQ.png"
    ),
    plot = p_AQ,
    width = 7,
    height = 5,
    dpi = 300
  )
}

