#' Spatial persistence among sampling events
#'
#' Calculates pairwise Spearman correlations among the spatial concentration
#' patterns observed during different sampling events, then returns the median
#' correlation for each event. The input must already be in site-by-event form.
#'
#' @param concentration_matrix Numeric matrix with sites in rows and sampling
#'   events in columns.
#' @param use Missing-value method passed to [stats::cor()].
#'
#' @return A named numeric vector with one spatial-persistence value per event.
#' @export
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 2, 4, 6), ncol = 2)
#' colnames(x) <- c("event_1", "event_2")
#' spatial_persistence(x)
spatial_persistence <- function(
    concentration_matrix,
    use = "pairwise.complete.obs"
) {
  validate_numeric_matrix(concentration_matrix, "concentration_matrix")
  if (ncol(concentration_matrix) < 2L) {
    stop("At least two sampling events are required.", call. = FALSE)
  }

  correlations <- suppressWarnings(stats::cor(
    concentration_matrix,
    method = "spearman",
    use = use
  ))
  diag(correlations) <- NA_real_

  apply(correlations, 1L, function(x) {
    if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  })
}

# Internal helper for calculate_spatial_persistence(): the Spearman
# correlation between two sampling events' spatial concentration patterns,
# restricted to sites sampled on both events.
calculate_event_correlation <- function(
    data,
    event_1,
    event_2,
    concentration_col,
    site_col,
    event_col,
    min_shared_sites = 5L
) {
  required_cols <- c(concentration_col, site_col, event_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  if (length(min_shared_sites) != 1L || is.na(min_shared_sites) ||
      min_shared_sites < 2 || min_shared_sites != as.integer(min_shared_sites)) {
    stop("min_shared_sites must be one integer of at least 2.", call. = FALSE)
  }

  keep <- data[[event_col]] %in% c(event_1, event_2)
  paired <- data[keep, required_cols, drop = FALSE]
  names(paired) <- c("concentration", "site", "event")
  paired <- paired[complete.cases(paired), , drop = FALSE]

  duplicate_keys <- duplicated(paired[c("site", "event")])
  if (any(duplicate_keys)) {
    stop("Duplicate site-event observations were found.", call. = FALSE)
  }

  first <- paired[paired$event == event_1, c("site", "concentration"), drop = FALSE]
  second <- paired[paired$event == event_2, c("site", "concentration"), drop = FALSE]
  names(first)[2] <- "concentration_event_1"
  names(second)[2] <- "concentration_event_2"
  shared <- merge(first, second, by = "site", all = FALSE, sort = FALSE)

  n_shared <- nrow(shared)
  enough_variation <- n_shared >= min_shared_sites &&
    length(unique(shared$concentration_event_1)) > 1L &&
    length(unique(shared$concentration_event_2)) > 1L

  rho <- if (enough_variation) {
    suppressWarnings(stats::cor(
      shared$concentration_event_1,
      shared$concentration_event_2,
      method = "spearman"
    ))
  } else {
    NA_real_
  }

  data.frame(
    event_1 = as.character(event_1),
    event_2 = as.character(event_2),
    n_shared_sites = n_shared,
    spearman_rho = rho,
    stringsAsFactors = FALSE
  )
}

#' Spatial persistence from a long data frame (legacy data-frame API)
#'
#' Calculates pairwise Spearman correlations among the spatial concentration
#' patterns observed during different sampling events, starting from a long
#' data frame rather than an already-built site-by-event matrix, then returns
#' the median correlation for each event alongside the full pairwise detail.
#'
#' @param data A long-format data frame with one row per site-event
#'   observation.
#' @param concentration_col Name of the concentration column.
#' @param site_col Name of the site column.
#' @param event_col Name of the sampling-event column.
#' @param min_shared_sites Minimum number of shared sites required to
#'   calculate a correlation for a pair of events; pairs with fewer shared
#'   sites (or without variation in either event) get `NA`.
#' @param Constituent Optional label (e.g. an analyte/constituent name)
#'   stamped onto every returned table as a `Constituent` column. Lets
#'   results from separate calls (one per analyte) be combined with
#'   `dplyr::bind_rows()` for a color-coded or faceted plot, the same way
#'   the hand-rolled `calc_SPpairs()` tags its own output. Default `NULL`
#'   leaves the tables unchanged (no column added).
#'
#' @return A list with three data frames: `by_event` (one row per event, with
#'   its median spatial persistence), `pairwise` (one row per event pair),
#'   and `event_comparisons` (the pairwise table duplicated in both
#'   directions, used to build `by_event`). Each carries a `Constituent`
#'   column when the `Constituent` argument is supplied.
#' @export
#'
#' @examples
#' d <- data.frame(
#'   event = rep(c("e1", "e2", "e3"), each = 3),
#'   site = rep(c("s1", "s2", "s3"), times = 3),
#'   conc = c(1, 2, 3, 2, 4, 5, 3, 5, 8)
#' )
#' calculate_spatial_persistence(d, "conc", "site", "event", min_shared_sites = 2)
calculate_spatial_persistence <- function(
    data,
    concentration_col,
    site_col,
    event_col,
    min_shared_sites = 5L,
    Constituent = NULL
) {
  column_args <- c(concentration_col, site_col, event_col)
  if (any(lengths(list(concentration_col, site_col, event_col)) != 1L) ||
      any(is.na(column_args)) || any(column_args == "")) {
    stop("Each column argument must be one non-missing column name.", call. = FALSE)
  }
  missing_cols <- setdiff(column_args, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(data[[concentration_col]])) {
    stop("The concentration column must be numeric.", call. = FALSE)
  }

  complete <- !is.na(data[[site_col]]) & !is.na(data[[event_col]]) &
    nzchar(trimws(as.character(data[[event_col]]))) &
    !is.na(data[[concentration_col]])
  clean <- data[complete, column_args, drop = FALSE]
  if (any(duplicated(clean[c(site_col, event_col)]))) {
    stop("Duplicate site-event observations were found.", call. = FALSE)
  }

  events <- unique(as.character(clean[[event_col]]))
  if (length(events) < 2L) {
    stop("At least two sampling events are required.", call. = FALSE)
  }

  pairs <- utils::combn(events, 2L, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    calculate_event_correlation(
      data = clean,
      event_1 = pair[1],
      event_2 = pair[2],
      concentration_col = concentration_col,
      site_col = site_col,
      event_col = event_col,
      min_shared_sites = min_shared_sites
    )
  }))
  rownames(pairwise) <- NULL

  comparisons <- rbind(
    data.frame(
      event = pairwise$event_1,
      comparison_event = pairwise$event_2,
      n_shared_sites = pairwise$n_shared_sites,
      spearman_rho = pairwise$spearman_rho,
      stringsAsFactors = FALSE
    ),
    data.frame(
      event = pairwise$event_2,
      comparison_event = pairwise$event_1,
      n_shared_sites = pairwise$n_shared_sites,
      spearman_rho = pairwise$spearman_rho,
      stringsAsFactors = FALSE
    )
  )

  by_event <- do.call(rbind, lapply(events, function(current_event) {
    event_rows <- comparisons[comparisons$event == current_event, , drop = FALSE]
    valid <- !is.na(event_rows$spearman_rho)
    data.frame(
      event = current_event,
      spatial_persistence = if (any(valid)) stats::median(event_rows$spearman_rho[valid]) else NA_real_,
      n_event_comparisons = sum(valid),
      total_possible_comparisons = nrow(event_rows),
      median_shared_sites = if (any(valid)) stats::median(event_rows$n_shared_sites[valid]) else NA_real_,
      min_shared_sites = if (any(valid)) min(event_rows$n_shared_sites[valid]) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(by_event) <- NULL

  if (!is.null(Constituent)) {
    pairwise$Constituent <- Constituent
    comparisons$Constituent <- Constituent
    by_event$Constituent <- Constituent
  }

  list(
    by_event = by_event,
    pairwise = pairwise,
    event_comparisons = comparisons
  )
}

#' Average spatial persistence across sampling events
#'
#' Calculates the mean of the per-event spatial-persistence values returned
#' by [calculate_spatial_persistence()], ignoring events with a missing
#' value.
#'
#' @param spatial_persistence_result The list returned by
#'   [calculate_spatial_persistence()].
#'
#' @return A single numeric value: the mean spatial persistence across all
#'   events with a non-missing value.
#' @export
#'
#' @examples
#' d <- data.frame(
#'   event = rep(c("e1", "e2", "e3"), each = 3),
#'   site = rep(c("s1", "s2", "s3"), times = 3),
#'   conc = c(1, 2, 3, 2, 4, 5, 3, 5, 8)
#' )
#' result <- calculate_spatial_persistence(d, "conc", "site", "event", min_shared_sites = 2)
#' average_spatial_persistence(result)
average_spatial_persistence <- function(spatial_persistence_result) {
  by_event <- spatial_persistence_result$by_event
  if (is.null(by_event) || is.null(by_event$spatial_persistence)) {
    stop(
      "spatial_persistence_result must be the list returned by calculate_spatial_persistence().",
      call. = FALSE
    )
  }
  mean(by_event$spatial_persistence, na.rm = TRUE)
}

