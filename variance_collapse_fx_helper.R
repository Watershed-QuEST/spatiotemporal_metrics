# Dependency-free collapse calculations; safe to source standalone.
# External observation preparation and PELT detection: scripts/variance_breakpoints.R.

#' Variance-collapse threshold
#'
#' Compares adjacent segment variances using externally detected breakpoints.
#' Supply the exact prepared observations used by the detector. This function
#' does not filter or sort observations and does not require changepoint.
#'
#' @param concentration Finite numeric vector of prepared concentrations.
#' @param area Positive finite contributing areas in nondecreasing order,
#'   aligned with `concentration`. Preserve the detector's order within ties.
#' @param breakpoints Required strictly increasing integer positions of the
#'   last observation of each internal segment; exclude the terminal row.
#'   Use `integer(0)` for completed detection with no breakpoint. Check the
#'   detector status first; skipped detection is not a no-breakpoint result.
#'
#' @return A list containing `threshold`, the ordered and scaled values in
#'   `data`, all adjacent comparisons in `changepoints`, `breakpoints`,
#'   `threshold_next_area`, `threshold_midpoint`, and `status`. The threshold
#'   is the lower bounding area of the first variance decrease. Status is
#'   `no_changepoint`, `variance_collapse_detected`, or `variance_increase_only`
#'   (no decrease, including equal variances). Each segment needs at least two
#'   observations. Keep the detector fit and settings in the source output.
#' @export
variance_collapse <- function(
    concentration,
    area,
    breakpoints
) {
  if (missing(breakpoints)) {
    stop("breakpoints is required; run breakpoint detection first.", call. = FALSE)
  }
  if (!is.numeric(concentration) || !is.numeric(area) ||
      !is.null(dim(concentration)) || !is.null(dim(area))) {
    stop("concentration and area must be numeric.", call. = FALSE)
  }
  if (length(concentration) != length(area)) {
    stop("concentration and area must have equal lengths.", call. = FALSE)
  }
  if (any(!is.finite(concentration)) || any(!is.finite(area)) ||
      any(area <= 0) || is.unsorted(area)) {
    stop("Supply finite prepared observations with positive areas in increasing order; do not filter or sort after detection.", call. = FALSE)
  }
  ordered <- data.frame(concentration = concentration, area = area)
  n <- nrow(ordered)
  if (!is.numeric(breakpoints) || !is.null(dim(breakpoints)) ||
      any(!is.finite(breakpoints)) || any(breakpoints != floor(breakpoints)) ||
      any(breakpoints < 1 | breakpoints >= n) || any(diff(breakpoints) <= 0)) {
    stop("breakpoints must be strictly increasing integer positions from 1 to n - 1; use integer(0) for no breakpoint.", call. = FALSE)
  }
  indices <- as.integer(breakpoints)
  if (any(diff(c(0L, indices, n)) < 2L)) {
    stop("Every segment must contain at least two observations.", call. = FALSE)
  }
  concentration_sd <- stats::sd(ordered$concentration)
  if (!is.finite(concentration_sd) || concentration_sd == 0) {
    stop("Concentration must have non-zero variance.", call. = FALSE)
  }
  ordered$scaled_concentration <- as.numeric(scale(ordered$concentration))

  ends <- c(indices, nrow(ordered))
  starts <- c(1L, indices + 1L)
  ordered$variance_segment <- rep(seq_along(ends), ends - starts + 1L)
  variances <- vapply(seq_along(ends), function(i) {
    stats::var(ordered$scaled_concentration[starts[i]:ends[i]])
  }, numeric(1))
  changes <- data.frame(
    index = indices,
    area = ordered$area[indices],
    next_area = ordered$area[indices + 1L],
    midpoint = (ordered$area[indices] + ordered$area[indices + 1L]) / 2,
    variance_before = variances[-length(variances)],
    variance_after = variances[-1L]
  )
  changes$variance_ratio <- changes$variance_after / changes$variance_before
  changes$is_variance_collapse <- changes$variance_after < changes$variance_before
  collapse <- which(changes$is_variance_collapse)

  list(
    threshold = if (length(collapse)) changes$area[collapse[1L]] else NA_real_,
    data = ordered,
    changepoints = changes,
    breakpoints = indices,
    threshold_next_area = if (length(collapse)) changes$next_area[collapse[1L]] else NA_real_,
    threshold_midpoint = if (length(collapse)) changes$midpoint[collapse[1L]] else NA_real_,
    status = if (!length(indices)) "no_changepoint" else if (length(collapse))
      "variance_collapse_detected" else "variance_increase_only"
  )
}

#' Variance collapse for grouped external detections
#'
#' Applies [variance_collapse()] to each completed detection from the script
#' helper `detect_variance_breakpoints_by_group()`. No detector is run here.
#' @param segmentation Grouped helper output with `groups` identifiers and
#'   aligned `results`, each containing prepared `data`, `breakpoints`, and status.
#' @return A list with `by_event` (one row per watershed/event/constituent) and
#'   `results` (aligned collapse lists, or NULL for skipped detections).
#'   Detector fit, settings, and source-row mappings remain in `segmentation`.
#' @export
variance_collapse_by_group <- function(segmentation) {
  if (!is.list(segmentation) || !is.data.frame(segmentation$groups) ||
      !is.list(segmentation$results) ||
      nrow(segmentation$groups) != length(segmentation$results) ||
      !all(c("watershed", "event", "constituent") %in% names(segmentation$groups))) {
    stop("Supply grouped output from detect_variance_breakpoints_by_group().", call. = FALSE)
  }
  results <- vector("list", length(segmentation$results))
  summaries <- vector("list", length(results))
  for (i in seq_along(results)) {
    detection <- segmentation$results[[i]]
    if (!is.list(detection) || !is.character(detection$status) ||
        length(detection$status) != 1L || is.na(detection$status) ||
        !detection$status %in% c("breakpoints_detected", "no_changepoint",
                                "too_few_sites", "zero_concentration_variance") ||
        !is.data.frame(detection$data) ||
        !all(c("concentration", "area") %in% names(detection$data)) ||
        !is.numeric(detection$breakpoints)) {
      stop("Invalid grouped detection at position ", i, ".", call. = FALSE)
    }
    completed <- detection$status %in% c("breakpoints_detected", "no_changepoint")
    if ((detection$status == "breakpoints_detected") != (length(detection$breakpoints) > 0L)) {
      stop("Detection status and breakpoints disagree at position ", i, ".", call. = FALSE)
    }
    collapse <- if (completed) variance_collapse(
      detection$data$concentration, detection$data$area, detection$breakpoints
    ) else NULL
    results[i] <- list(collapse)
    summaries[[i]] <- data.frame(
      n_sites = nrow(detection$data), detection_status = detection$status,
      status = if (completed) collapse$status else detection$status,
      collapse_detected = completed && collapse$status == "variance_collapse_detected",
      threshold = if (completed) collapse$threshold else NA_real_,
      threshold_next_area = if (completed) collapse$threshold_next_area else NA_real_,
      threshold_midpoint = if (completed) collapse$threshold_midpoint else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  summary <- if (length(summaries)) do.call(rbind, summaries) else data.frame(
    n_sites = integer(), detection_status = character(), status = character(),
    collapse_detected = logical(), threshold = numeric(),
    threshold_next_area = numeric(), threshold_midpoint = numeric()
  )
  list(by_event = cbind(segmentation$groups, summary), results = results)
}

