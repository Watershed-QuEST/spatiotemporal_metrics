#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-08-03

# This script fits the real synoptic data provided in 01_build_spatial_data.R, validates the fit, and produces everything the metric-specific scripts need for their sensitivity analyses.
# The statistical model producing synthetic data (used throughout) is: log(concentration) = SSN2's kriged spatial mean + a deterministic seasonal shift + spatially-correlated noise (optionally also correlated across campaigns; see Part C's "red noise"). See spatiotemporal_helpers.R for the shared math.

# This script is organized as:
#   PART A -- clean the toy dataset, assign campaign IDs to consecutive-day campaigns
#   PART B -- fit SSN2 spatial covariance models (tail-up/tail-down/Euclidean/nugget) to each constituent's campaign-averaged data, with diagnostics (Torgegrams, AICc + leave-one-out cross-validation (LOOCV) comparison) and semi-automated model selection
#   PART C -- fit a log-normal area + seasonal (sin/cos DOY) model per constituent, and check its residuals for leftover temporal autocorrelation once SSN2's spatial mean and the seasonal cycle are both removed. This informs the red-noise parameter in Part E.
#   PART D -- save a reusable set up per analyte/constituent that saves SSN2's kriged mean + spatial Cholesky factor, the seasonal coefficients, and a temporal autocorrelation estimate so that downstream scripts can each draw fresh Monte Carlo fields without recomputing the computationally-intensive Cholesky step. See notes in this section for an explanation of the Cholesky factor!
#   PART E -- generate and save ONE static extended synthetic dataset. We want this to use for metrics that need more data to demonstrate than the real toy dataset alone. 

# Requires: 01_build_spatial_data.R must be run first, and spatiotemporal_helpers.R must be saved in the same folder as this script.

# Outputs:
#   data/nm_clean.csv: toy data cleaned up of non-campaign data
#   data/nm_field_setups.rds: per-constituent set up needed for downstream scripts to draw fresh MC fields from
#   data/nm_synthetic_extended.csv: the static synthetic extended dataset (long format)
#   data/*.csv, plots/*.png: diagnostics for each section

#### Packages ####
library(tidyverse)
library(sf)
library(SSN2)

# raster package messes with dplyr::select functionality in some sessions;
# this keeps select() resolved to dplyr's version regardless.
select <- dplyr::select

# connects this script to the helper functions in spatiotemporal_helpers.R
source("spatiotemporal_helpers.R")


#### Read in data ####

# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download files
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM-BR Toy dataset.csv"],
                            path = "drivedata/toy.csv",
                            overwrite = T)
googledrive::drive_download(file = toy_files$id[toy_files$name=="All sites.xlsx"], 
                            path = "drivedata/All sites.xlsx",
                            overwrite = T)

#### Configure in/outputs and file structure ####

data_dir  <- "."
data_file <- "drivedata/toy.csv"
out_dir   <- "."

plot_dir     <- file.path(out_dir, "plots")
data_out_dir <- file.path(out_dir, "data")
for (d in c(plot_dir, data_out_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## This whole pipeline (scripts 2-4) is scoped to one watershed at a time. To adapt to a different single watershed, change these and the matching entries in 01_build_spatial_data.R's `networks` list.
project          <- "nm"
ssn_path         <- "SantaFe.ssn"
ssn_site_id_col  <- "Site"      # matches the toy dataset directly
ssn_additive_col <- "afvFlow"   # real flow-accumulation AFV, from 01's Part B2
excluded_sites   <- c("USF01", "USF02", "USF30")  # must match 01_build_spatial_data.R
constituents     <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

max_gap_days   <- 1   # max day-gap allowed within one campaign run
min_run_length <- 2   # runs shorter than this are discarded as one-off revisits

n_synthetic_months <- 60      # 5 years of monthly campaigns for the static extended dataset (Part E)
time_unit_days     <- 30.44   # days per synthetic "month"; also the time unit temporal_rho below is expressed per

set.seed(42)

#### PART A -- Clean toy dataset ####

# Note: nm's sites can't all be visited in a day, so a campaign is a run of consecutive calendar dates (gap <= max_gap_days). Isolated dates not part of a run of >= min_run_length are one-off revisits, not full campaigns, and are discarded.

raw <- read_csv(file.path(data_dir, data_file), show_col_types = FALSE) %>%
  mutate(Date = mdy(Date), DOY = yday(Date)) %>%
  filter(Project == project, !Site %in% excluded_sites)

date_runs <- raw %>%
  distinct(Date) %>%
  arrange(Date) %>%
  mutate(gap = as.numeric(Date - lag(Date)),
         new_run = is.na(gap) | gap > max_gap_days,
         run_id = cumsum(new_run)) %>%
  add_count(run_id, name = "run_length")

kept_dates <- date_runs %>%
  filter(run_length >= min_run_length) %>%
  mutate(CampaignID = paste0(project, "_C", dense_rank(run_id))) %>%
  select(Date, CampaignID)

clean <- raw %>%
  inner_join(kept_dates, by = "Date") %>%   # drops one-off dates automatically
  arrange(CampaignID, Site) %>%
  mutate(logArea = log(Area.m2))

# examine data before saving!
View(clean)
nrow(raw)
nrow(clean)
print(clean %>% count(CampaignID))

write_csv(clean, file.path(data_out_dir, "nm_clean.csv"))

#### PART B -- Fit SSN2 spatial stream network covariance model ####

# Read in ssn object
ssn_object <- readRDS(file.path(data_out_dir, "ssn_objects.rds"))[[project]]

## ssn_lm() reads its response/predictor columns directly from ssn.object$obs, so the cleaned data needs to be pivoted (one column per CampaignID per constituent) and joined onto that sf object.
wide <- clean %>%
  select(Site, logArea, CampaignID, all_of(constituents)) %>%
  pivot_wider(id_cols = c(Site, logArea), names_from = CampaignID,
              values_from = all_of(constituents), names_glue = "{.value}__{CampaignID}")
ssn_object$obs <- ssn_object$obs %>% left_join(wide, by = setNames("Site", ssn_site_id_col))

## This produces diagnostic torgegrams (semivariance by flow-connected, flow-unconnected, and Euclidean distance) per campaign. 
# NOTE: Plots provide a starting point for intuition about covariance structure, for reference. Note that actual model selection occurs below via AICc + LOOCV comparison. By using model comparison over torgrgrams as diagnostic, we are taking the stance here that synthetic data and MC analyses should incorporate stream network spatial structuring even if the torgegrams don't scream that there is a strong or specific spatial structure in the data. By taking this stance, we can generate synthetic data with realistic spatial variance, which is important given that we are interested in characterizing spatial variance! So, take a look at these, but we don't recommend using them to decide what model structure to go with alone.
camps <- unique(clean$CampaignID)
n_torgegrams <- 0
for (cc in constituents) {
  for (camp in camps) {
    resp <- paste(cc, camp, sep = "__")
    if (!resp %in% names(ssn_object$obs)) next
    if (sum(!is.na(ssn_object$obs[[resp]])) < 2) next  # Torgegram needs >= 2 non-NA sites

    form <- as.formula(paste0("log(`", resp, "`) ~ logArea"))
    tg <- tryCatch(
      Torgegram(formula = form, ssn.object = ssn_object, type = c("flowcon", "flowuncon", "euclid")),
      error = function(e) { message("Torgegram failed for ", resp, ": ", e$message); NULL }
    )
    if (!is.null(tg)) {
      png(file.path(plot_dir, paste0("Torgegram_", resp, ".png")), width = 900, height = 600)
      plot(tg); dev.off()
      n_torgegrams <- n_torgegrams + 1
    }
  }
}

## This defines a list of candidate covariance structures for model comparison. 
##   none          -- no spatial structure at all: nugget only (AICc baseline)
##   tailup        -- tail-up exponential only (stream-network mixing)
##   euclid_only   -- Euclidean exponential only (landscape-scale/elevation driver, no stream-network mixing)
##   tailup_euclid -- tail-up exponential + Euclidean exponential
covariance_specs <- list(
  none          = list(tailup_type = "none",       taildown_type = "none", euclid_type = "none"),
  tailup        = list(tailup_type = "exponential", taildown_type = "none", euclid_type = "none"),
  euclid_only   = list(tailup_type = "none",       taildown_type = "none", euclid_type = "exponential"),
  tailup_euclid = list(tailup_type = "exponential", taildown_type = "none", euclid_type = "exponential")
)

## This helper function sidesteps a bug in ssn_lm() where fitting any tailup covariance against a response column containing NA throws "non-conformable matrix dimensions in additive * (b == 0)". none/ euclid_only tolerate NA responses fine (they drop those rows internally, same as base lm() would); tailup's additive-weight ("branching") apparently isn't subsetted to match when rows get dropped, causing the dimension mismatch. Fix: use ssn_subset() to properly drop NA-response sites Returns the original ssn.object unchanged if there's nothing to drop.
drop_na_response <- function(ssn.object, resp, label) {
  n_na <- sum(is.na(ssn.object$obs[[resp]]))
  if (n_na == 0) return(ssn.object)

  ssn.object <- ssn_update_path(ssn.object, normalizePath(ssn.object$path, mustWork = FALSE))

  ## Scratch object, written to tempdir() (not the Dropbox-synced project
  ## folder) -- disposable input to ssn_lm(), not a deliverable.
  sub_path <- file.path(tempdir(), paste0(
    tools::file_path_sans_ext(basename(ssn_path)), "_", make.names(label), "_sub.ssn"
  ))
  subset_expr <- bquote(!is.na(.(as.symbol(resp))))

  sub_ssn <- tryCatch(
    do.call(ssn_subset, list(ssn = ssn.object, path = sub_path,
                              subset = subset_expr, clip = FALSE, overwrite = TRUE)),
    error = function(e) { message(resp, ": ssn_subset failed: ", e$message); NULL }
  )
  if (is.null(sub_ssn)) return(NULL)

  ok <- tryCatch({
    ssn_create_distmat(sub_ssn, predpts = "predpts", among_predpts = TRUE, overwrite = TRUE)
    TRUE
  }, error = function(e) { message(resp, ": ssn_create_distmat on subset failed: ", e$message); FALSE })
  if (!ok) return(NULL)

  cat("  ", resp, "-- dropped", n_na, "NA-response site(s) before fitting\n")
  sub_ssn
}

## Fit the covariance_specs to each constituent's mean concentration across its available campaigns. Averaging can sharpen a spatial signal that looks like noise in any single campaign. This uses the arithmetic mean of raw concentration per site.
all_fits    <- list()
all_site_ids <- list()
varcomp_rows <- list()

# Fit the models + calculate AICc and leave-one-out cross-validation (LOOCV)
for (cc in constituents) {
  camp_cols <- intersect(paste(cc, camps, sep = "__"), names(ssn_object$obs))
  if (length(camp_cols) == 0) next

  resp_avg <- paste(cc, "avg", sep = "__")
  obs_vals <- st_drop_geometry(ssn_object$obs)[, camp_cols, drop = FALSE]
  avg_vals <- rowMeans(as.matrix(obs_vals), na.rm = TRUE)
  avg_vals[is.nan(avg_vals)] <- NA
  ssn_object$obs[[resp_avg]] <- avg_vals

  fit_ssn <- drop_na_response(ssn_object, resp_avg, resp_avg)
  if (is.null(fit_ssn)) next

  form <- as.formula(paste0("log(`", resp_avg, "`) ~ logArea"))

  for (model_name in names(covariance_specs)) {
    spec <- covariance_specs[[model_name]]

    ## do.call() (not a direct ssn_lm(...) call) is deliberate
    fit <- tryCatch(
      do.call(ssn_lm, c(list(formula = form, ssn.object = fit_ssn, additive = ssn_additive_col), spec)),
      error = function(e) { message(resp_avg, " [", model_name, "] failed: ", e$message); NULL }
    )
    if (is.null(fit)) next

    key <- paste(resp_avg, model_name, sep = "__")
    all_fits[[key]] <- fit
    ## Site IDs in the exact row order ssn_lm() saw them -- needed in Part C to attach loocv()'s per-row predictions back to a real Site.
    all_site_ids[[key]] <- st_drop_geometry(fit_ssn$obs)[[ssn_site_id_col]]

    ## This adds in leave-one-out cross-validation (LOOCV), in addition to AICc. LOOCV's provides a direct check of whether the extra covariance flexibility actually buys better predictions, not just a parsimony evaluation. 
    cv <- tryCatch(loocv(fit), error = function(e) {
      message(resp_avg, " [", model_name, "] loocv failed: ", e$message); NULL
    })

    vc <- varcomp(fit) %>%
      mutate(Constituent = cc, Model = model_name, AICc = AICc(fit), Key = key, .before = 1)
    if (!is.null(cv)) vc <- bind_cols(vc, cv[rep(1, nrow(vc)), , drop = FALSE])
    varcomp_rows[[key]] <- vc

    cat(sprintf("  %-30s %-13s AICc=%7.2f", resp_avg, model_name, AICc(fit)),
        if (!is.null(cv)) sprintf("  RMSPE=%.3f cor2=%.3f cover.95=%.2f", cv$RMSPE, cv$cor2, cv$cover.95)
        else "  (loocv failed)", "\n")
  }
}

varcomp_summary <- bind_rows(varcomp_rows)
write_csv(varcomp_summary, file.path(data_out_dir, "SSN_covariance_summary.csv"))
saveRDS(all_fits, file.path(data_out_dir, "SSN_fitted_models.rds"))

## SEMI-AUTOMATED MODEL SELECTION: combine AICc + LOOCV into one rule
# AICc is a valid basis for comparing these directly: the response and fixed-effect formula (log(conc) ~ logArea) are identical across all four, only the covariance structure differs.
# Selection rule, applied per constituent:
#   1. Does ANY spatial structure beat "no structure" by a meaningful AICc margin? If not,  "none" (no spatial structure) is a legitimate, simulatable result in its own right. We recommend adding it to covariance_choice_overrides below to run it through next steps.
#   2. Among specs that clear step 1, drop any whose LOOCV calibration looks broken (std.MSPE far from 1, or cover.95 well below nominal). This indicates overconfident uncertainty estimates.
#   3. Among what's left, lowest RMSPE (tie-break: highest cor2) wins, UNLESS it's within rmspe_tie_tol of another spec. User must break a real tie using the covariance_choice_overrides tibble below.
# Note: Thresholds below are judgment calls -- loosen/tighten here if clearly too strict/permissive.
select_best_covariance_fit <- function(varcomp_summary,
                                        aicc_window = 2,
                                        min_cover95 = 0.70,
                                        std_mspe_bounds = c(0.3, 3),
                                        rmspe_tie_tol = 0.02,
                                        overrides = NULL) {
  if (is.null(overrides)) overrides <- tibble(Constituent = character(), Model = character())

  by_fit <- varcomp_summary %>%
    distinct(Constituent, Model, Key, AICc, bias, std.bias, MSPE, RMSPE,
             std.MSPE, RAV, cor2, cover.80, cover.90, cover.95)

  none_ref <- by_fit %>% filter(Model == "none") %>%
    group_by(Constituent) %>% summarize(none_aicc = min(AICc), .groups = "drop")

  by_fit %>%
    filter(Model != "none") %>%
    left_join(none_ref, by = "Constituent") %>%
    group_by(Constituent) %>%
    group_modify(function(df, key) {
      none_aicc <- df$none_aicc[1]
      manual <- overrides %>% filter(Constituent == key$Constituent)

      ## A manual override wins
      if (nrow(manual) == 1) {
        if (manual$Model[1] == "none") {
          none_row <- by_fit %>% filter(Constituent == key$Constituent, Model == "none") %>%
            slice_min(AICc, n = 1, with_ties = FALSE) %>%
            mutate(none_aicc = none_aicc) %>% select(-Constituent)
          if (nrow(none_row) == 1) {
            cat("  ", key$Constituent, "-- using manual override: 'none' (no spatial structure).\n")
            return(none_row)
          }
          cat("  ", key$Constituent, "-- override requests 'none' but no 'none' fit was found -- ignoring it.\n")
        } else {
          manual_row <- df %>% filter(Model == manual$Model[1])
          if (nrow(manual_row) == 1) {
            cat("  ", key$Constituent, "-- using manual override:", manual$Model[1],
                "(bypassing the automated AICc/LOOCV selection rule below).\n")
            return(manual_row)
          }
          cat("  ", key$Constituent, "-- override requests '", manual$Model[1],
              "' but no such fit was found for this constituent -- ignoring it.\n", sep = "")
        }
      }

      best_aicc <- min(df$AICc)
      if (is.na(none_aicc) || best_aicc > none_aicc - aicc_window) {
        cat("  ", key$Constituent, "-- no spatial structure beats 'none' by >=", aicc_window,
            "AICc (best spatial AICc =", round(best_aicc, 2), ", none AICc =", round(none_aicc, 2),
            ") -- 'none' is a legitimate result here, not just a disqualifier. Add a row to",
            "covariance_choice_overrides to run it through anyway. Skipping for now.\n")
        return(df[0, ])
      }

      calibrated <- df %>% filter(std.MSPE >= std_mspe_bounds[1], std.MSPE <= std_mspe_bounds[2],
                                   cover.95 >= min_cover95)
      pool <- if (nrow(calibrated) > 0) calibrated else df
      if (nrow(calibrated) == 0) {
        cat("  ", key$Constituent, "-- no spec passed the calibration check -- falling back to",
            "best-RMSPE anyway; check its LOOCV numbers before trusting it.\n")
      }

      pool <- pool %>% arrange(RMSPE, desc(cor2))
      tied <- pool %>% filter(RMSPE <= pool$RMSPE[1] * (1 + rmspe_tie_tol))
      if (nrow(tied) == 1) return(tied)

      cat("\n  ", key$Constituent, "-- TIE:", nrow(tied), "spec(s) within",
          paste0(round(rmspe_tie_tol * 100), "%"), "RMSPE of each other -- not auto-resolving:\n")
      print(tied %>% select(Model, AICc, RMSPE, cor2, std.MSPE, cover.95))
      cat("  Add a row to covariance_choice_overrides to pick one. Skipping for now.\n\n")
      df[0, ]
    }) %>%
    ungroup() %>%
    rename(key = Key)
}

# Two-call workflow, by design: the first call has no overrides
best_fits <- select_best_covariance_fit(varcomp_summary)

# ENTER MANUAL MODEL CHOICES HERE. They will be reviewed against the printout above.
covariance_choice_overrides <- tibble::tribble(
  ~Constituent,    ~Model,
  "NPOC..mg.C.L.", "tailup_euclid",
  "TDN..mg.N.L.",  "tailup"
)

# Call again with overrides
best_fits <- select_best_covariance_fit(varcomp_summary, overrides = covariance_choice_overrides)

# Review the covariance model each constituent will be simulated from:
print(best_fits %>% select(Constituent, Model, AICc, RMSPE, cor2))

#### PART C -- Fit seasonally-varying model and check residual temporal autocorrelation ####
# The basic model is:
# log(C) = b0 + b1*log(Area) + b2*sin(2*pi*DOY/365) + b3*cos(2*pi*DOY/365) + e

seasonal_models <- list()
for (cc in constituents) {
  sub <- clean[!is.na(clean[[cc]]) & clean[[cc]] > 0, ]
  fit_df <- data.frame(y = log(sub[[cc]]), logArea = sub$logArea,
                        sin_doy = sin(2 * pi * sub$DOY / 365), cos_doy = cos(2 * pi * sub$DOY / 365))
  fit <- lm(y ~ logArea + sin_doy + cos_doy, data = fit_df)
  seasonal_models[[cc]] <- list(fit = fit, beta = unname(coef(fit)), sigma = summary(fit)$sigma,
                                 area_rng = range(sub$logArea), doy_rng = range(sub$DOY), obs = sub)
}

# Fitted area+seasonal model coefficients (Intercept, logArea, sin, cos) and sigma (e):
for (cc in constituents) {
  m <- seasonal_models[[cc]]
  cat(sprintf("  %-15s beta=[%6.3f %6.3f %6.3f %6.3f]  sigma=%.3f\n",
              cc, m$beta[1], m$beta[2], m$beta[3], m$beta[4], m$sigma))
}

## Fit-check plots -- check that the deterministic area+season curve tracks the real data at all 
# NOTE: these estimates are before spatial structure from SSN2 is incorporated, so you should expect a lot of spread that can be attributed to spatial heterogeneity. The goal is not a perfect fit, but to see if any systematice pattern over time in the real data is being represented in predictions. Loops over analytes/constituents and saves plot to file. 
for (cc in constituents) {
  m <- seasonal_models[[cc]]
  area_seq <- seq(m$area_rng[1], m$area_rng[2], length.out = 50)
  area_curve <- data.frame(logArea = area_seq, sin_doy = 0, cos_doy = 0) %>%
    mutate(pred = predict(m$fit, newdata = ., se.fit = TRUE)$fit,
           se = predict(m$fit, newdata = ., se.fit = TRUE)$se.fit,
           Conc = exp(pred), Conc_lo = exp(pred - 1.96 * se), Conc_hi = exp(pred + 1.96 * se))
  p <- ggplot(m$obs, aes(x = exp(logArea), y = .data[[cc]])) +
    geom_point(alpha = 0.5) +
    geom_ribbon(data = area_curve, aes(x = exp(logArea), y = NULL, ymin = Conc_lo, ymax = Conc_hi), alpha = 0.2) +
    geom_line(data = area_curve, aes(x = exp(logArea), y = Conc)) +
    labs(x = "Drainage area (m2)", y = cc, title = paste(cc, "-- fit check: area"),
         subtitle = "Line + band = fitted mean and 95% CI at seasonally-neutral DOY; points = observed") +
    theme_bw()
  ggsave(file.path(plot_dir, paste0("Fit_check_area_", cc, ".png")), p, width = 8, height = 6, dpi = 150)

  doy_seq <- seq(m$doy_rng[1], m$doy_rng[2], length.out = 100)
  season_curve <- data.frame(logArea = mean(m$obs$logArea), sin_doy = sin(2 * pi * doy_seq / 365),
                              cos_doy = cos(2 * pi * doy_seq / 365), DOY = doy_seq) %>%
    mutate(pred = predict(m$fit, newdata = ., se.fit = TRUE)$fit,
           se = predict(m$fit, newdata = ., se.fit = TRUE)$se.fit,
           Conc = exp(pred), Conc_lo = exp(pred - 1.96 * se), Conc_hi = exp(pred + 1.96 * se))
  p <- ggplot(m$obs, aes(x = DOY, y = .data[[cc]])) +
    geom_point(alpha = 0.5) +
    geom_ribbon(data = season_curve, aes(x = DOY, y = NULL, ymin = Conc_lo, ymax = Conc_hi), alpha = 0.2) +
    geom_line(data = season_curve, aes(x = DOY, y = Conc)) +
    labs(x = "Day of year", y = cc, title = paste(cc, "-- fit check: season"),
         subtitle = "Line + band = fitted mean and 95% CI at mean observed log(Area); points = observed") +
    theme_bw()
  ggsave(file.path(plot_dir, paste0("Fit_check_seasonal_", cc, ".png")), p, width = 8, height = 6, dpi = 150)
}

## Residual temporal autocorrelation calculation (temporal_rho) and diagnostic
# This checks the assumption Part D/E's temporal autocorrelation parameter (temporal_rho) is built on: once variation attributable to SSN2's spatial mean and the seasonal shift are both removed, is what's left of a site's residual still correlated from one campaign to the next? Ideally no, but below we pass on the remaining temporal autocorrelation as a temporal_rho parameter in the model so that synthetic data reproduces it and thus shows realistic temporal structure. 
# This uses SSN2's own leave-one-out prediction (so a site can't "predict itself") to represent spatial patterns so that they don't show up as temporal autocorrelation. The workflow pairs up temporally adjacent residuals within each site so we can look for temporal autocorrelation across more pairs than any single site's 1-3 campaigns could support alone. The resulting temporal_rho is a low-power estimate of the true temporal autocorrelation given that only 3 campaigns exist, but still provides a reasonable estimate for synthesizing data with realistic temporal structures when seasonal patterns do not adequately represent temporal structure. 

resid_lag1_pairs <- pmap_dfr(
  list(best_fits$Constituent, best_fits$Model, best_fits$key),
  function(cc, model_name, key) {
    fit <- all_fits[[key]]
    site_ids <- all_site_ids[[key]]
    if (is.null(fit) || is.null(site_ids)) return(tibble())

    cv <- tryCatch(loocv(fit, cv_predict = TRUE), error = function(e) {
      cat("  ", cc, ": loocv(cv_predict = TRUE) failed --", e$message, "\n"); NULL
    })
    if (is.null(cv) || is.null(cv$cv_predict)) return(tibble())
    mu_site_loocv <- setNames(cv$cv_predict, site_ids)

    beta <- seasonal_models[[cc]]$beta
    ref_doy <- mean(clean$DOY, na.rm = TRUE)

    obs <- seasonal_models[[cc]]$obs
    obs$mu_loocv <- mu_site_loocv[obs$Site]
    obs <- obs[!is.na(obs$mu_loocv), ]  
    if (nrow(obs) == 0) return(tibble())

    obs$Residual <- log(obs[[cc]]) - obs$mu_loocv - seasonal_shift(obs$DOY, beta, ref_doy)

    obs %>%
      arrange(Site, Date) %>%
      group_by(Site) %>%
      filter(n() >= 2) %>%
      mutate(Residual_prev = lag(Residual)) %>%
      ungroup() %>%
      filter(!is.na(Residual_prev)) %>%
      transmute(Constituent = cc, Site, Date, Residual, Residual_prev)
  })

resid_autocorr_summary <- tibble(Constituent = constituents, n_pairs = 0L, lag1_cor = NA_real_)
if (nrow(resid_lag1_pairs) >= 4) {
  resid_autocorr_summary <- resid_lag1_pairs %>%
    group_by(Constituent) %>%
    summarize(n_pairs = n(), lag1_cor = cor(Residual, Residual_prev), .groups = "drop")

  cat("  Pooled lag-1 (consecutive-campaign, same-site) residual correlation:\n")
  print(resid_autocorr_summary)

  write_csv(resid_lag1_pairs, file.path(data_out_dir, "Residual_lag1_pairs.csv"))
  write_csv(resid_autocorr_summary, file.path(data_out_dir, "Residual_lag1_summary.csv"))

  p_resid <- ggplot(resid_lag1_pairs, aes(x = Residual_prev, y = Residual)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    facet_wrap(~ Constituent, scales = "free") +
    labs(x = "Residual, previous campaign (same site)", y = "Residual, this campaign",
         title = "Residual temporal autocorrelation check") +
    theme_bw()
  p_resid
  ggsave(file.path(plot_dir, "Residual_lag1_plot.png"), p_resid, width = 8, height = 6, dpi = 150)
} else {
  cat("  Not enough consecutive same-site campaign pairs (with a usable LOOCV mean) to check.\n")
}

#### PART D -- Save set up for simulating data ####

# see "spatiotemporal_helpers.R for details about numerical stability if you get a warning message about the build_field_setup() function

# Note: the Cholesky factor used here is essentially a "square root" of a covariance matrix: for a matrix X describing how a set of variables should correlate, there's a matrix L such that multiplying L by its own transpose gives back X. This is useful here because multiplying independent noise by L (the spatial covariance) is the trick that "mixes" random noise with spatially-correlated noise. It's easy to generate plain, independent random noise (rnorm()), but not easy to directly generate noise that's already correlated in some specific pattern - the Cholesky factor addresses this issue. 

temporal_rho_overrides <- list()  # e.g. list(`NPOC..mg.C.L.` = 0.2) to override a specific constituent

field_setups <- list()
for (i in seq_len(nrow(best_fits))) {
  cc <- best_fits$Constituent[i]
  key <- best_fits$key[i]
  fit <- all_fits[[key]]
  if (is.null(fit)) next

  setup <- build_field_setup(fit)

  rho <- temporal_rho_overrides[[cc]]
  if (is.null(rho)) {
    est <- resid_autocorr_summary$lag1_cor[resid_autocorr_summary$Constituent == cc]
    ## draw_field()'s exponential-decay temporal covariance can only represent rho between 0 and 1, never negative correlation; rho = 0 (white noise in time).
    rho <- if (length(est) == 1 && !is.na(est)) max(0, min(0.9, est)) else 0
  }

  field_setups[[cc]] <- list(setup = setup, beta = seasonal_models[[cc]]$beta,
                              ref_doy = mean(clean$DOY, na.rm = TRUE), temporal_rho = rho,
                              model = best_fits$Model[i])
  cat("  ", cc, "-- field setup saved (m =", setup$m, "predpts, model =", best_fits$Model[i],
      ", temporal_rho =", round(rho, 3), ")\n")
}

saveRDS(field_setups, file.path(data_out_dir, "nm_field_setups.rds"))

#### PART E -- Generate and save static extended synthetic dataset ####
# this does one realization (not a Monte Carlo) covering every predicted location x n_synthetic_months monthly campaigns, for use for metrics that need more data to demonstrate. Actual sensitivity tests redraw fresh fields from field_setups above instead of resampling this file, so their uncertainty reflects real variable uncertainty. 

t_days <- seq(15, by = time_unit_days, length.out = n_synthetic_months)  # ~mid-month
doy_seq <- t_days %% 365.25

synthetic_rows <- list()
for (cc in names(field_setups)) {
  fs <- field_setups[[cc]]
  noise <- draw_field(fs$setup, t_days, rho = fs$temporal_rho, time_unit_days = time_unit_days)
  shift_mat <- matrix(seasonal_shift(doy_seq, fs$beta, fs$ref_doy),
                       nrow = fs$setup$m, ncol = n_synthetic_months, byrow = TRUE)
  field <- exp(fs$setup$mu_cond + shift_mat + noise)

  colnames(field) <- paste0("SynthCamp", seq_len(n_synthetic_months))
  synthetic_rows[[cc]] <- as_tibble(field) %>%
    mutate(Constituent = cc, predID = paste0("P", seq_len(fs$setup$m)), .before = 1) %>%
    pivot_longer(starts_with("SynthCamp"), names_to = "CampaignNum", values_to = "Conc") %>%
    mutate(CampaignNum = as.integer(str_remove(CampaignNum, "SynthCamp")),
           Day = t_days[CampaignNum], DOY = doy_seq[CampaignNum])
}

synthetic_extended <- bind_rows(synthetic_rows)
write_csv(synthetic_extended, file.path(data_out_dir, "nm_synthetic_extended.csv"))

## Check plot: the synthetic concentrations themselves (mean across all predpts locations each month, with a light band showing the 25th-75th percentile spread across locations that month), over the full 5 years -- eyeballs whether the seasonal cycle and red-noise addition (if temporal_rho > 0) look like a plausible concentration record.
p_synth_check <- synthetic_extended %>%
  group_by(Constituent, CampaignNum, Day) %>%
  summarize(Mean_Conc = mean(Conc), P25_Conc = quantile(Conc, 0.25),
            P75_Conc = quantile(Conc, 0.75), .groups = "drop") %>%
  ggplot(aes(x = Day / 365.25, y = Mean_Conc)) +
  geom_ribbon(aes(ymin = P25_Conc, ymax = P75_Conc), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 1) +
  labs(x = "Synthetic year", y = "Concentration (mean across predpts, that month)",
       title = "Static extended dataset")+
  theme_bw()
ggsave(file.path(plot_dir, "Synthetic_extended_mean_by_month.png"), p_synth_check, width = 8, height = 6, dpi = 150)


