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
library(ggplot2)
library(patchwork)


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

#### Establish Independence Classification ####

# This is used for sensitivity analyses later on
# Classification was established by knowledge of the sites and their hydrologic connectivity.

# Independent NM sites
nm_independent <- c(
  "USF04", "USF13", "USF16", "USF17", "USF18",
  "USF15", "USF22", "USF23", "USF24", "USF25",
  "USF26", "USF28", "USF29"
)

toy <- toy %>%
  mutate(
    Project = str_to_upper(Project),
    Independence = case_when(
      
      # BR classification
      Project == "BR" &
        (Site == "BRM07" | !str_detect(str_to_upper(Site), "M")) ~
        "Independent",
      
      Project == "BR" ~ "Non-independent",
      
      # NM classification
      Project == "NM" & Site %in% nm_independent ~
        "Independent",
      
      Project == "NM" ~ "Non-independent",
      
      TRUE ~ NA_character_
    )
  )

#===============================================================================
# SITE-LEVEL SUMMARY
#===============================================================================

# Calculate the mean of log-transformed mass loads and discharge
# for each site. These site-level means are used to fit the regressions.

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

# This creates 8 separate models, 4 for each catchment, 2 for each solute (NPOC and TDN), 
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
# PLOTS: SAR POINT PLOTS BY WATERSHED AND SOLUTE
#===============================================================================

# Y-axis limits
ylims <- list(
  NM = c(-1.5, 0.8),
  BR = c(-1.5, 0.8)
)

# Solutes to plot
solutes <- c("NPOC", "TDN")

# Store plots
plots_A <- list()
plots_AQ <- list()

for (solute in solutes) {
  
  # Determine column prefix
  solute_prefix <- paste0("log", solute)
  
  for (w in watersheds) {
    
    #---------------------------------------------------------------------------
    # COLUMN NAMES
    #---------------------------------------------------------------------------
    
    # Area-only SARs
    sar_A_mean <- paste0(
      "SAR_", solute_prefix, "_A_", w, "_mean"
    )
    
    sar_A_se <- paste0(
      "SAR_", solute_prefix, "_A_", w, "_se"
    )
    
    # Area + discharge SARs
    sar_AQ_mean <- paste0(
      "SAR_", solute_prefix, "_AQ_", w, "_mean"
    )
    
    sar_AQ_se <- paste0(
      "SAR_", solute_prefix, "_AQ_", w, "_se"
    )
    
    # Filter watershed data
    plot_data <- residual_summary %>%
      filter(Project == w)
    
    
    #---------------------------------------------------------------------------
    # AREA-ONLY MODEL
    #---------------------------------------------------------------------------
    
    p_A <- ggplot(plot_data) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "black"
      ) +
      
      geom_point(
        aes(
          x = reorder(Site, .data[[sar_A_mean]]),
          y = .data[[sar_A_mean]]
        ),
        color = "skyblue4",
        size = 3
      ) +
      
      geom_errorbar(
        aes(
          x = reorder(Site, .data[[sar_A_mean]]),
          ymin = .data[[sar_A_mean]] - .data[[sar_A_se]],
          ymax = .data[[sar_A_mean]] + .data[[sar_A_se]]
        ),
        width = 0.2,
        colour = "orange",
        alpha = 0.9,
        linewidth = 1.1
      ) +
      
      coord_cartesian(
        ylim = ylims[[w]]
      ) +
      
      labs(
        x = "Site",
        y = expression("Mean SAR (" * log[10] * " scale)"),
        title = paste0(w, " ", solute, " SAR (A Model)")
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
    
    plots_A[[paste(solute, w, sep = "_")]] <- p_A
    
    
    #---------------------------------------------------------------------------
    # AREA + DISCHARGE MODEL
    #---------------------------------------------------------------------------
    
    p_AQ <- ggplot(plot_data) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "black"
      ) +
      
      geom_point(
        aes(
          x = reorder(Site, .data[[sar_AQ_mean]]),
          y = .data[[sar_AQ_mean]]
        ),
        color = "skyblue4",
        size = 3
      ) +
      
      geom_errorbar(
        aes(
          x = reorder(Site, .data[[sar_AQ_mean]]),
          ymin = .data[[sar_AQ_mean]] - .data[[sar_AQ_se]],
          ymax = .data[[sar_AQ_mean]] + .data[[sar_AQ_se]]
        ),
        width = 0.2,
        colour = "orange",
        alpha = 0.9,
        linewidth = 1.1
      ) +
      
      coord_cartesian(
        ylim = ylims[[w]]
      ) +
      
      labs(
        x = "Site",
        y = expression("Mean SAR (" * log[10] * " scale)"),
        title = paste0(w, " ", solute, " SAR (AQ Model)")
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
    
    plots_AQ[[paste(solute, w, sep = "_")]] <- p_AQ
    
  }
}


#===============================================================================
# COMBINE WATERSHEDS SIDE BY SIDE AND SAVE
#===============================================================================

for (solute in solutes) {
  
  #---------------------------------------------------------------------------
  # AREA-ONLY MODEL: NM + BR
  #---------------------------------------------------------------------------
  
  p_A_combined <-
    plots_A[[paste(solute, "NM", sep = "_")]] +
    plots_A[[paste(solute, "BR", sep = "_")]] +
    plot_layout(ncol = 2)
  
  print(p_A_combined)
  
  ggsave(
    filename = paste0("plots/SAR_", solute, "_A.png"),
    plot = p_A_combined,
    width = 12,
    height = 5,
    dpi = 300
  )
  
  
  #---------------------------------------------------------------------------
  # AREA + DISCHARGE MODEL: NM + BR
  #---------------------------------------------------------------------------
  
  p_AQ_combined <-
    plots_AQ[[paste(solute, "NM", sep = "_")]] +
    plots_AQ[[paste(solute, "BR", sep = "_")]] +
    plot_layout(ncol = 2)
  
  print(p_AQ_combined)
  
  ggsave(
    filename = paste0("plots/SAR_", solute, "_AQ.png"),
    plot = p_AQ_combined,
    width = 12,
    height = 5,
    dpi = 300
  )
  
}


#===============================================================================
# MODEL COMPARISON
#===============================================================================

get_metrics <- function(model) {
  
  model_summary <- summary(model)
  
  tibble(
    R_squared = model_summary$r.squared,
    Adjusted_R_squared = model_summary$adj.r.squared,
    AIC = AIC(model),
    RMSE = sqrt(mean(residuals(model)^2)),
    n = nobs(model)
  )
}

model_comparison <- list()

for (w in watersheds) {
  
  for (s in c("logNPOC", "logTDN")) {
    
    model_A <- lm_A[[w]][[s]]
    model_AQ <- lm_AQ[[w]][[s]]
    
    metrics_A <- get_metrics(model_A)
    metrics_AQ <- get_metrics(model_AQ)
    
    model_comparison[[paste(w, s, sep = "_")]] <- tibble(
      
      Project = w,
      Solute = s,
      
      n = metrics_A$n,
      
      A_R_squared = metrics_A$R_squared,
      AQ_R_squared = metrics_AQ$R_squared,
      
      A_adj_R_squared = metrics_A$Adjusted_R_squared,
      AQ_adj_R_squared = metrics_AQ$Adjusted_R_squared,
      
      A_AIC = metrics_A$AIC,
      AQ_AIC = metrics_AQ$AIC,
      
      A_RMSE = metrics_A$RMSE,
      AQ_RMSE = metrics_AQ$RMSE,
      
      Delta_AIC = metrics_AQ$AIC - metrics_A$AIC,
      
      Delta_adj_R_squared =
        metrics_AQ$Adjusted_R_squared -
        metrics_A$Adjusted_R_squared,
      
      Q_F_test_p = anova(model_A, model_AQ)$`Pr(>F)`[2]
    )
  }
}

model_comparison <- bind_rows(model_comparison)

print(model_comparison, width = Inf)


#===============================================================================
# SENSITIVITY TESTS
#===============================================================================

#===============================================================================
# SENSITIVITY TEST 1: BACK-TRANSFORM SAR RESIDUALS
#===============================================================================

# Back-transform observation-level SAR residuals (leads to a mathematical ratio)
# SAR ratios represent observed mass load / predicted mass load

sar_cols <- grep(
  "^SAR_",
  names(toy_log),
  value = TRUE
)

toy_log <- toy_log %>%
  mutate(
    across(
      all_of(sar_cols),
      ~ 10^.x,
      .names = "Ratio_{.col}"
    )
  )

#===============================================================================
# SUMMARIZE BACK-TRANSFORMED SAR RATIOS
#===============================================================================

# Standardized ratio column names (without watershed suffix)

ratio_names <- c(
  "Ratio_SAR_logNPOC_A",
  "Ratio_SAR_logTDN_A",
  "Ratio_SAR_logNPOC_AQ",
  "Ratio_SAR_logTDN_AQ"
)


# Combine NM and BR values into shared columns

for (nm in ratio_names) {
  
  toy_log[[nm]] <- ifelse(
    toy_log$Project == "NM",
    toy_log[[paste0(nm, "_NM")]],
    toy_log[[paste0(nm, "_BR")]]
  )
  
}


# Functions for calculating means and standard errors

mean_or_na <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}

se_or_na <- function(x) {
  n <- sum(!is.na(x))
  
  if (n < 2) {
    NA_real_
  } else {
    sd(x, na.rm = TRUE) / sqrt(n)
  }
}


# Summarize ratios by site

ratio_summary <- toy_log %>%
  group_by(Site, Project, Area.m2) %>%
  summarise(
    across(
      all_of(ratio_names),
      list(
        mean = mean_or_na,
        se = se_or_na
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

View(ratio_summary)

#XXXXXXX

#===============================================================================
# SENSITIVITY ANALYSIS: BACK-TRANSFORMED SAR RATIOS
#===============================================================================

solutes <- c("NPOC", "TDN")

plots_ratio_A <- list()
plots_ratio_AQ <- list()

for (solute in solutes) {
  
  solute_prefix <- paste0("log", solute)
  
  for (w in watersheds) {
    
    # Column names for back-transformed ratios
    ratio_A_mean <- paste0("Ratio_SAR_", solute_prefix, "_A_mean")
    ratio_A_se   <- paste0("Ratio_SAR_", solute_prefix, "_A_se")
    
    ratio_AQ_mean <- paste0("Ratio_SAR_", solute_prefix, "_AQ_mean")
    ratio_AQ_se   <- paste0("Ratio_SAR_", solute_prefix, "_AQ_se")
    
    # Filter to watershed
    plot_data <- ratio_summary %>%
      filter(Project == w)
    
    #-------------------------------
    # A MODEL
    #-------------------------------
    
    p_A <- ggplot(plot_data) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "black"
      ) +
      geom_point(
        aes(
          x = reorder(Site, .data[[ratio_A_mean]]),
          y = .data[[ratio_A_mean]]
        ),
        color = "skyblue4",
        size = 3
      ) +
      geom_errorbar(
        aes(
          x = reorder(Site, .data[[ratio_A_mean]]),
          ymin = pmax(0, .data[[ratio_A_mean]] - .data[[ratio_A_se]]),
          ymax = .data[[ratio_A_mean]] + .data[[ratio_A_se]]
        ),
        width = 0.2,
        color = "orange",
        alpha = 0.9,
        linewidth = 1.1
      ) +
      labs(
        x = "Site",
        y = "Mean observed / predicted mass-load ratio",
        title = paste0(w, " ", solute, " SAR Ratio (A Model)")
      ) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    plots_ratio_A[[paste(solute, w, sep = "_")]] <- p_A
    
    
    #-------------------------------
    # AQ MODEL
    #-------------------------------
    
    p_AQ <- ggplot(plot_data) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "black"
      ) +
      geom_point(
        aes(
          x = reorder(Site, .data[[ratio_AQ_mean]]),
          y = .data[[ratio_AQ_mean]]
        ),
        color = "skyblue4",
        size = 3
      ) +
      geom_errorbar(
        aes(
          x = reorder(Site, .data[[ratio_AQ_mean]]),
          ymin = pmax(0, .data[[ratio_AQ_mean]] - .data[[ratio_AQ_se]]),
          ymax = .data[[ratio_AQ_mean]] + .data[[ratio_AQ_se]]
        ),
        width = 0.2,
        color = "orange",
        alpha = 0.9,
        linewidth = 1.1
      ) +
      labs(
        x = "Site",
        y = "Mean observed / predicted mass-load ratio",
        title = paste0(w, " ", solute, " SAR Ratio (AQ Model)")
      ) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    plots_ratio_AQ[[paste(solute, w, sep = "_")]] <- p_AQ
  }
}


#===============================================================================
# COMBINE AND SAVE SENSITIVITY ANALYSIS PLOTS
#===============================================================================

for (solute in solutes) {
  
  # A model: NM + BR
  p_A_combined <-
    plots_ratio_A[[paste(solute, "NM", sep = "_")]] +
    plots_ratio_A[[paste(solute, "BR", sep = "_")]] +
    plot_layout(ncol = 2)
  
  print(p_A_combined)
  
  ggsave(
    filename = paste0(
      "plots/SAR_log_sensitivity_", solute, "_A.png"
    ),
    plot = p_A_combined,
    width = 12,
    height = 5,
    dpi = 300
  )
  
  
  # AQ model: NM + BR
  p_AQ_combined <-
    plots_ratio_AQ[[paste(solute, "NM", sep = "_")]] +
    plots_ratio_AQ[[paste(solute, "BR", sep = "_")]] +
    plot_layout(ncol = 2)
  
  print(p_AQ_combined)
  
  ggsave(
    filename = paste0(
      "plots/SAR_log_sensitivity_", solute, "_AQ.png"
    ),
    plot = p_AQ_combined,
    width = 12,
    height = 5,
    dpi = 300
  )
}

# Quick interpretation, values close to 1 indicated exporting as expected based on regressions created,
# 0-1 indicates exporting less than expected, and >1 indicates exporting more than expected.
print(p_A_combined)
print(p_AQ_combined)


#===============================================================================
# SENSITIVITY TEST 2: AREA VS. DISCHARGE VS. AREA + DISCHARGE
#===============================================================================

#-------------------------------------------------------------------------------
# 1. CALCULATE Q-ONLY SAR RESIDUALS
#-------------------------------------------------------------------------------

solutes <- c("NPOC", "TDN")

for (w in watersheds) {
  
  for (solute in solutes) {
    
    mass_col <- paste0("log", solute)
    sar_col  <- paste0("SAR_", mass_col, "_Q_", w)
    
    # Select observations from the watershed
    idx <- which(toy_log$Project == w)
    
    # Calculate site-averaged log values for model fitting
    site_means <- toy_log[idx, ] %>%
      group_by(Site) %>%
      summarise(
        logQ = mean(logQ, na.rm = TRUE),
        log_mass = mean(.data[[mass_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(
        is.finite(logQ),
        is.finite(log_mass)
      )
    
    # Fit Q-only model
    q_model <- lm(log_mass ~ logQ, data = site_means)
    
    # Predict mass load for individual observations
    predicted_mass <- predict(
      q_model,
      newdata = data.frame(logQ = toy_log$logQ[idx])
    )
    
    # Calculate observation-level Q-only SAR
    toy_log[[sar_col]] <- NA_real_
    
    toy_log[[sar_col]][idx] <-
      toy_log[[mass_col]][idx] - predicted_mass
  }
}


#-------------------------------------------------------------------------------
# 2. COMPARE SITE-LEVEL MEAN SAR VALUES
#-------------------------------------------------------------------------------

sar_model_comparison <- purrr::map_dfr(watersheds, function(w) {
  
  purrr::map_dfr(solutes, function(solute) {
    
    a_col  <- paste0("SAR_log", solute, "_A_", w)
    q_col  <- paste0("SAR_log", solute, "_Q_", w)
    aq_col <- paste0("SAR_log", solute, "_AQ_", w)
    
    toy_log %>%
      filter(Project == w) %>%
      group_by(Project, Site) %>%
      summarise(
        A_SAR = mean(.data[[a_col]], na.rm = TRUE),
        Q_SAR = mean(.data[[q_col]], na.rm = TRUE),
        AQ_SAR = mean(.data[[aq_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Solute = solute)
    
  })
  
})

# View site-level comparison
View(sar_model_comparison)


#-------------------------------------------------------------------------------
# 3. CREATE PAIRWISE MODEL COMPARISONS
#-------------------------------------------------------------------------------

sar_pairwise <- bind_rows(
  
  sar_model_comparison %>%
    transmute(
      Project, Site, Solute,
      Comparison = "A vs. Q",
      First_SAR = A_SAR,
      Second_SAR = Q_SAR
    ),
  
  sar_model_comparison %>%
    transmute(
      Project, Site, Solute,
      Comparison = "A vs. AQ",
      First_SAR = A_SAR,
      Second_SAR = AQ_SAR
    ),
  
  sar_model_comparison %>%
    transmute(
      Project, Site, Solute,
      Comparison = "Q vs. AQ",
      First_SAR = Q_SAR,
      Second_SAR = AQ_SAR
    )
)


#-------------------------------------------------------------------------------
# 4. CORRELATION AND AGREEMENT BETWEEN MODELS
#-------------------------------------------------------------------------------

sar_correlation <- sar_pairwise %>%
  group_by(Project, Solute, Comparison) %>%
  summarise(
    
    n = sum(complete.cases(First_SAR, Second_SAR)),
    
    Pearson_r = cor(
      First_SAR,
      Second_SAR,
      method = "pearson",
      use = "complete.obs"
    ),
    
    Spearman_rho = cor(
      First_SAR,
      Second_SAR,
      method = "spearman",
      use = "complete.obs"
    ),
    
    Sign_agreement_pct = mean(
      sign(First_SAR) == sign(Second_SAR),
      na.rm = TRUE
    ) * 100,
    
    .groups = "drop"
  )

print(sar_correlation, width = Inf)


#-------------------------------------------------------------------------------
# 5. PLOT PAIRWISE SAR COMPARISONS
#-------------------------------------------------------------------------------

p_sar_comparison <- ggplot(
  sar_pairwise,
  aes(x = First_SAR, y = Second_SAR)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(
    color = "skyblue4",
    size = 3,
    alpha = 0.85
  ) +
  facet_grid(Project ~ Solute + Comparison) +
  labs(
    x = "SAR from first model",
    y = "SAR from second model",
    title = "Sensitivity of SAR to Model Predictors"
  ) +
  theme_classic(base_size = 13)

print(p_sar_comparison)

ggsave(
  "plots/SAR_A_Q_sensitivity.png",
  plot = p_sar_comparison,
  width = 14,
  height = 8,
  dpi = 300
)




#===============================================================================
# SENSITIVITY TEST 1: INDEPENDENT VS. NON-INDEPENDENT SITES
#===============================================================================

# Only calculate the two sensitivity groups.
# The original all-sites SAR values will come from residual_summary.

site_sets <- c("Independent", "Non-independent")
solutes <- c("NPOC", "TDN")
models <- c("A", "AQ")


#-------------------------------------------------------------------------------
# 1. Calculate SAR for independent and non-independent sites
#-------------------------------------------------------------------------------

calculate_independence_sar <- function(
    data,
    site_set,
    watershed,
    solute,
    model_type
) {
  
  mass_col <- paste0("log", solute)
  
  # Filter to watershed and site independence group
  data_sub <- data %>%
    filter(
      Project == watershed,
      Independence == site_set
    )
  
  # Select model predictors
  predictors <- if (model_type == "A") {
    "logArea"
  } else {
    c("logArea", "logQ")
  }
  
  # Calculate site-level means for model fitting
  site_means <- data_sub %>%
    group_by(Site) %>%
    summarise(
      across(
        all_of(c(predictors, mass_col)),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    filter(
      if_all(
        all_of(c(predictors, mass_col)),
        is.finite
      )
    )
  
  # Check that there are enough sites to fit the model
  if (nrow(site_means) <= length(predictors)) {
    warning(
      paste(
        "Insufficient sites:",
        site_set,
        watershed,
        solute,
        model_type
      )
    )
    
    return(tibble())
  }
  
  # Fit SAR regression
  model_formula <- reformulate(
    predictors,
    response = mass_col
  )
  
  sar_model <- lm(
    model_formula,
    data = site_means
  )
  
  # Calculate observation-level SAR
  data_sub %>%
    mutate(
      Predicted_log_mass = predict(
        sar_model,
        newdata = .
      ),
      Observed_log_mass = .data[[mass_col]],
      SAR = Observed_log_mass - Predicted_log_mass,
      Site_Set = site_set,
      Model = model_type,
      Solute = solute
    ) %>%
    select(
      Project,
      Site_Set,
      Independence,
      Site,
      Solute,
      Model,
      Observed_log_mass,
      Predicted_log_mass,
      SAR
    )
}


# Run models for both watersheds, solutes, and model types
independence_sar_obs <- tidyr::crossing(
  Site_Set = site_sets,
  Project = watersheds,
  Solute = solutes,
  Model = models
) %>%
  purrr::pmap_dfr(
    function(Site_Set, Project, Solute, Model) {
      
      calculate_independence_sar(
        data = toy_log,
        site_set = Site_Set,
        watershed = Project,
        solute = Solute,
        model_type = Model
      )
    }
  )


#-------------------------------------------------------------------------------
# 2. Summarize SAR by site
#-------------------------------------------------------------------------------

independence_sar_summary <- independence_sar_obs %>%
  group_by(
    Project,
    Site_Set,
    Independence,
    Site,
    Solute,
    Model
  ) %>%
  summarise(
    Mean_SAR = mean(SAR, na.rm = TRUE),
    SE_SAR = if (sum(is.finite(SAR)) > 1) {
      sd(SAR, na.rm = TRUE) / sqrt(sum(is.finite(SAR)))
    } else {
      NA_real_
    },
    n = sum(is.finite(SAR)),
    .groups = "drop"
  )


# View site-level sensitivity results
View(independence_sar_summary)


#-------------------------------------------------------------------------------
# 3. Summarize independent and non-independent groups
#-------------------------------------------------------------------------------

independence_model_summary <- independence_sar_summary %>%
  group_by(Project, Site_Set, Solute, Model) %>%
  summarise(
    Number_of_sites = n_distinct(Site),
    Mean_SAR = mean(Mean_SAR, na.rm = TRUE),
    SD_SAR = sd(Mean_SAR, na.rm = TRUE),
    Median_SAR = median(Mean_SAR, na.rm = TRUE),
    .groups = "drop"
  )

print(independence_model_summary, n = Inf, width = Inf)


#-------------------------------------------------------------------------------
# 4. Extract original all-sites SAR values
#-------------------------------------------------------------------------------

# The original SAR calculation is already stored in residual_summary.
# Extract the site-level mean SAR columns and convert them to long format.

original_sar_long <- residual_summary %>%
  select(Project, Site, matches("^SAR_log(NPOC|TDN)_(A|AQ)_(NM|BR)_mean$")) %>%
  pivot_longer(
    cols = matches("^SAR_log(NPOC|TDN)_(A|AQ)_(NM|BR)_mean$"),
    names_to = c("Solute", "Model", "Project_from_name"),
    names_pattern = "^SAR_log(NPOC|TDN)_(A|AQ)_(NM|BR)_mean$",
    values_to = "Mean_SAR"
  ) %>%
  mutate(
    Site_Set = "All sites",
    Project = as.character(Project)
  ) %>%
  filter(
    Project == Project_from_name
  ) %>%
  select(
    Project,
    Site,
    Site_Set,
    Solute,
    Model,
    Mean_SAR
  )


#-------------------------------------------------------------------------------
# 5. Combine original and sensitivity SAR values
#-------------------------------------------------------------------------------

# Keep the same structure for all three groups.

independence_boxplot_data <- bind_rows(
  
  # Original all-sites SAR values
  original_sar_long,
  
  # Independent and non-independent SAR values
  independence_sar_summary %>%
    select(
      Project,
      Site,
      Site_Set,
      Solute,
      Model,
      Mean_SAR
    )
  
) %>%
  mutate(
    Site_Set = factor(
      Site_Set,
      levels = c(
        "Independent",
        "Non-independent",
        "All sites"
      )
    ),
    Project = factor(Project, levels = c("NM", "BR")),
    Solute = factor(Solute, levels = c("NPOC", "TDN")),
    Model = factor(Model, levels = c("A", "AQ"))
  ) %>%
  filter(is.finite(Mean_SAR))


# Check number of sites in each group
independence_boxplot_data %>%
  group_by(Project, Site_Set, Solute, Model) %>%
  summarise(
    Number_of_sites = n_distinct(Site),
    .groups = "drop"
  ) %>%
  print(n = Inf)


#-------------------------------------------------------------------------------
# 6. Create boxplots
#-------------------------------------------------------------------------------

p_independence <- ggplot(
  independence_boxplot_data,
  aes(
    x = Site_Set,
    y = Mean_SAR
  )
) +
  
  # Reference line at zero
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  
  # Boxplots
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    fill = "grey85",
    color = "black"
  ) +
  
  # Individual site values
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.7
  ) +
  
  # Separate panels by watershed, solute, and model
  facet_grid(Project ~ Solute + Model) +
  
  labs(
    title = "SAR Sensitivity to Site Independence",
    x = NULL,
    y = "Mean Site SAR"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )


# Display plot
print(p_independence)


#-------------------------------------------------------------------------------
# 7. Save plot
#-------------------------------------------------------------------------------

dir.create(
  "plots/",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "plots/SAR_independence_sensivity.png",
  plot = p_independence,
  width = 12,
  height = 8,
  dpi = 300
)
