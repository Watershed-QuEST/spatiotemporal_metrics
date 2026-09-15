# Script helpers: source this file; these are not package exports.
# See variance_breakpoints.md for the observation and breakpoint contract.

prepare_variance_observations <- function(concentration, area) {
  if (!is.numeric(concentration) || !is.numeric(area) ||
      !is.null(dim(concentration)) || !is.null(dim(area)) ||
      length(concentration) != length(area)) {
    stop("concentration and area must be aligned numeric vectors.", call. = FALSE)
  }
  keep <- is.finite(concentration) & is.finite(area)
  if (any(area[keep] <= 0)) {
    stop("Finite contributing areas must be positive.", call. = FALSE)
  }
  prepared <- data.frame(
    observation_index = which(keep),
    concentration = concentration[keep],
    area = area[keep]
  )
  prepared <- prepared[order(prepared$area, prepared$observation_index), , drop = FALSE]
  rownames(prepared) <- NULL
  concentration_sd <- stats::sd(prepared$concentration)
  prepared$scaled_concentration <- rep(NA_real_, nrow(prepared))
  if (is.finite(concentration_sd) && concentration_sd > 0) {
    prepared$scaled_concentration <-
      (prepared$concentration - mean(prepared$concentration)) / concentration_sd
  }
  prepared
}

detect_variance_breakpoints <- function(
    concentration, area, penalty = "MBIC", pen_value = 0,
    min_segment_length = 2L, min_sites = 6L
) {
  valid_integer <- function(x, minimum) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= minimum && x == floor(x) && x <= .Machine$integer.max
  }
  if (!valid_integer(min_segment_length, 2L)) {
    stop("min_segment_length must be one integer of at least 2.", call. = FALSE)
  }
  if (!valid_integer(min_sites, 4L)) {
    stop("min_sites must be one integer of at least 4.", call. = FALSE)
  }
  prepared <- prepare_variance_observations(concentration, area)
  n <- nrow(prepared)
  settings <- list(
    method = "PELT", test_stat = "Normal", penalty = penalty,
    pen_value = pen_value, min_segment_length = min_segment_length,
    min_sites = min_sites, minimum_required = max(min_sites, 2 * min_segment_length),
    standardization = "overall sample mean and sample SD",
    changepoint_version = NA_character_
  )
  result <- list(
    data = prepared, breakpoints = integer(0), fit = NULL,
    settings = settings, status = "too_few_sites",
    n_input = length(concentration), n_used = n,
    excluded_indices = setdiff(seq_along(concentration), prepared$observation_index),
    tied_areas = anyDuplicated(prepared$area) > 0L
  )
  if (n < settings$minimum_required) return(result)
  if (all(is.na(prepared$scaled_concentration))) {
    result$status <- "zero_concentration_variance"
    return(result)
  }
  if (!requireNamespace("changepoint", quietly = TRUE)) {
    stop("Install 'changepoint' to run breakpoint detection.", call. = FALSE)
  }
  fit <- changepoint::cpt.var(
    data = prepared$scaled_concentration, method = "PELT",
    test.stat = "Normal", penalty = penalty, pen.value = pen_value,
    minseglen = min_segment_length, class = TRUE, param.estimates = TRUE
  )
  indices <- changepoint::cpts(fit)
  indices <- sort(unique(as.integer(indices[indices > 0 & indices < n])))
  if (any(diff(c(0L, indices, n)) < min_segment_length)) {
    stop("Detector returned segments shorter than min_segment_length.", call. = FALSE)
  }
  result$breakpoints <- indices
  result$fit <- fit
  result$settings$changepoint_version <- as.character(utils::packageVersion("changepoint"))
  result$status <- if (length(indices)) "breakpoints_detected" else "no_changepoint"
  result
}

# Grouped, wide-table entry point, following the CV column-argument convention.
# event/watershed = NULL means the supplied table represents a single group
# on that dimension. Missing grouping labels are rejected rather than pooled.
detect_variance_breakpoints_by_group <- function(
    data, concentration, area, event = NULL, watershed = NULL,
    penalty = "MBIC", pen_value = 0, min_segment_length = 2L, min_sites = 6L
) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("data must be a nonempty data frame.", call. = FALSE)
  }
  column_name <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if (!is.character(concentration) || !length(concentration) ||
      anyNA(concentration) || any(!nzchar(concentration)) || !column_name(area) ||
      (!is.null(event) && !column_name(event)) ||
      (!is.null(watershed) && !column_name(watershed))) {
    stop("Supply concentration column names and single area/event/watershed column names.", call. = FALSE)
  }
  required <- c(unname(concentration), area, event, watershed)
  if (length(setdiff(required, names(data)))) {
    stop("Missing columns: ", paste(setdiff(required, names(data)), collapse = ", "), call. = FALSE)
  }
  labels <- names(concentration)
  if (is.null(labels)) labels <- unname(concentration)
  blank <- is.na(labels) | !nzchar(labels)
  labels[blank] <- unname(concentration)[blank]
  if (anyDuplicated(labels)) stop("Constituent labels must be unique.", call. = FALSE)
  keys <- data.frame(
    watershed = if (is.null(watershed)) rep("all", nrow(data)) else as.character(data[[watershed]]),
    event = if (is.null(event)) rep("all", nrow(data)) else as.character(data[[event]]),
    stringsAsFactors = FALSE
  )
  if (anyNA(keys) || any(!nzchar(trimws(keys$watershed))) ||
      any(!nzchar(trimws(keys$event)))) {
    stop("Resolve missing or blank watershed/event identifiers before detection.", call. = FALSE)
  }
  groups <- unique(keys)
  results <- list()
  identifiers <- list()
  for (i in seq_len(nrow(groups))) {
    rows <- which(keys$watershed == groups$watershed[i] & keys$event == groups$event[i])
    for (j in seq_along(concentration)) {
      detection <- detect_variance_breakpoints(
        data[[concentration[j]]][rows], data[[area]][rows],
        penalty = penalty, pen_value = pen_value,
        min_segment_length = min_segment_length, min_sites = min_sites
      )
      # Map local prepared positions back to rows of the full source table.
      detection$data$source_row <- rows[detection$data$observation_index]
      detection$excluded_source_rows <- rows[detection$excluded_indices]
      k <- length(results) + 1L
      results[[k]] <- detection
      identifiers[[k]] <- data.frame(
        watershed = groups$watershed[i], event = groups$event[i],
        constituent = labels[j], stringsAsFactors = FALSE
      )
    }
  }
  prepared_rows <- Map(function(id, result) {
    cbind(id[rep(1L, nrow(result$data)), , drop = FALSE], result$data)
  }, identifiers, results)
  structure(list(groups = do.call(rbind, identifiers), results = results,
                 data = do.call(rbind, prepared_rows)),
            class = c("variance_breakpoint_groups", "list"))
}
