#### READ ME

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-08-03

# This script contains a small set of shared functions needed to build synthetic data and perform Mote-Carlo sensitivity analyses. These are the pieces that are easy to mess up if the scripts are adapted, so are best kept in this separate script. 
#
# Requires: SSN2 (for covmatrix()/predict() inside build_field_setup()).

## seasonal_shift() - adds a seasonal cycle to data over time relative to a reference day-of-year:
# Part of the model used to synthesize data: log(C) = spatial mean + seasonal shift + noise. beta is the coefficient vector from an lm(y ~ logArea + sin_doy + cos_doy) fit; only the sin/cos terms (beta[3], beta[4]) are used here. The spatial terms (intercept, logArea) are deliberately NOT reused, since SSN2's own fitted mean (mu_cond, from build_field_setup() below) already covers the spatial part of the model.

seasonal_shift <- function(doy, beta, ref_doy) {
  term <- function(d) beta[3] * sin(2 * pi * d / 365) + beta[4] * cos(2 * pi * d / 365)
  term(doy) - term(ref_doy)
}

## build_field_setup() - conditional mean for a fitted SSN2 model, evaluated at "predpts"
# Standard kriging/Gaussian-process result: given fitted covariance matrices Sigma_oo (obs-obs), Sigma_po (pred-obs), and Sigma_pp (pred-pred), the kriged (conditional) mean at the prediction points is mu_cond (what predict() returns), and the conditional covariance of the predictions given the observed data is Sigma_cond = Sigma_pp - Sigma_po %*% solve(Sigma_oo) %*% t(Sigma_po).  # L_space below is chol(Sigma_cond) (R's convention: upper-triangular, t(L_space) %*% L_space = Sigma_cond)
# draw_field() uses t(L_space) as the matrix that turns iid noise into spatially correlated noise with covariance Sigma_cond.
# Note that this is a computationally intensive step (called once per analyte) to build synthetic data, then saved to disk, so later scripts can each just readRDS() it and draw fresh fields without repeating.

build_field_setup <- function(fit) {
  mu_cond  <- as.numeric(predict(fit, newdata = "predpts"))
  Sigma_oo <- covmatrix(fit)
  Sigma_po <- covmatrix(fit, newdata = "predpts", cov_type = "pred.obs")
  Sigma_pp <- covmatrix(fit, newdata = "predpts", cov_type = "pred.pred")

  Sigma_cond <- Sigma_pp - Sigma_po %*% solve(Sigma_oo, t(Sigma_po))
  m <- nrow(Sigma_cond)

  # "Ridge" for numerical stability - A covariance matrix like Sigma_cond should always be mathematically valid for chol(), but floating-point rounding error from preceding matrix math can leave it just barely invalid, causing chol() to fail (chol() is the thing that we draw correlated random values from). The fix is a "ridge": adding a tiny bit of independent variance to each point's own diagonal entry, which nudges the matrix back to valid without meaningfully changing what it represents. We scale that nudge to the matrix's own variance (rather than a fixed number) and retries with a 100x-bigger nudge, up to 6 times, if it still fails. This produces a warning so the user can proceed with caution, since a bigger-than-default nudge being needed is a warning sign that the underlying covariance fit is numerically fragile.
  ridge_scale <- mean(diag(Sigma_cond))
  ridge <- 1e-8 * ridge_scale
  L_space <- tryCatch(chol(Sigma_cond + diag(ridge, m)), error = function(e) NULL)
  attempt <- 1
  while (is.null(L_space) && attempt <= 6) {
    ridge <- ridge * 100
    attempt <- attempt + 1
    L_space <- tryCatch(chol(Sigma_cond + diag(ridge, m)), error = function(e) NULL)
  }
  if (is.null(L_space)) {
    stop("build_field_setup(): Sigma_cond is not positive definite even after inflating ",
         "the ridge to ", signif(ridge, 3), " (", attempt, " attempts). This usually means ",
         "the fitted covariance parameters are poorly identified -- if this model was forced ",
         "in via covariance_choice_overrides, that's a real warning sign about trusting it ",
         "for simulation, not just a numerical inconvenience. Inspect its tailup/euclid ",
         "range and nugget estimates (summary(fit)) before proceeding.")
  }
  if (ridge > 1e-8 * ridge_scale) {
    message("build_field_setup(): needed a larger-than-default ridge (", signif(ridge, 3),
            " vs. the usual ", signif(1e-8 * ridge_scale, 3), ") to make Sigma_cond positive ",
            "definite -- this fit's covariance may be numerically fragile; treat its ",
            "simulated points with some caution.")
  }

  list(mu_cond = mu_cond, L_space = L_space, m = m)
}

## draw_field() - does one fresh draw of spatially (and optionally temporally) correlated zero-mean noise, on a log scale to match log scale of predictive model. Returns an m x length(t_days) matrix. Each column (a single synthetic campaign) has spatial covariance Sigma_cond. We add setup$mu_cond and seasonal_shift() to this to get a full synthetic log-concentration field; draw_field() only returns the noise term, so the actual predictive model (mu_cond + shift + noise) stays visible in the main script instead of being hidden in here.
# rho is the correlation between two campaigns time_unit_days apart. If rho > 0, this additionally correlates DIFFERENT columns (campaigns) in time, via an exponential decay in elapsed time. rho = 0 (the default) reproduces independent-across-campaigns ("white") noise. When rho > 0 is added, it adds a lag-1 autocorrelation term (red noise) of the magnitude seen in the real data. Note that there is a limitation in this method that does not allow negative rhos because we are not providing enough data to estimate oscillations in time. Rather, as rho is set up here, we can represent how fast correlation decays from 1 toward 0 over time, which is a good approximation for producing realistic temporal noise into the synthesized data, but not exactly predicting it. 
# A reasonable starting point for rho, per constituent, is the pooled lag-1 residual correlation from the residual serial-correlation diagnostic in 02_build_synthetic_data.R (Residual_lag1_summary.csv), floored at 0.

draw_field <- function(setup, t_days, rho = 0, time_unit_days = 30.44) {
  if (rho < 0 || rho >= 1) {
    stop("draw_field(): rho must be in [0, 1) -- got ", rho, ". A negative value can't be ",
         "represented by this exponential-decay temporal covariance (see the comment above ",
         "this function); floor it at 0 before calling draw_field() instead of passing it through.")
  }

  n_t <- length(t_days)
  Z <- matrix(rnorm(setup$m * n_t), nrow = setup$m, ncol = n_t)
  A_space <- t(setup$L_space)     # A_space %*% t(A_space) = Sigma_cond
  spatial_noise <- A_space %*% Z  # m x n_t, spatially correlated, still
                                   # independent column-to-column (white in time)

  if (rho == 0 || n_t == 1) return(spatial_noise)

  dt <- abs(outer(t_days, t_days, "-")) / time_unit_days
  Sigma_time <- rho ^ dt
  L_time <- chol(Sigma_time + diag(1e-8, n_t))
  B_time <- t(L_time)             # B_time %*% t(B_time) = Sigma_time

  spatial_noise %*% t(B_time)
}
