# Standalone and package-ready CV helpers
#
# This file deliberately uses only base R and stats. It can be:
#   1. loaded as part of the watershedmetrics package, or
#   2. copied into another repository and loaded with source().


#' Spatial coefficient of variation
#'
#' Calculates the coefficient of variation across sites separately for each
#' sampling event using sample SD (denominator `n - 1`), then calculates the
#' unweighted mean of the eligible event CVs for each watershed.
#'
#' @param data A data frame containing concentration and event columns.
#' @param concentration Character string naming one numeric concentration
#'   column, or a named character vector mapping constituent labels to numeric
#'   concentration columns.
#' @param event Character string naming the sampling-event column.
#' @param watershed Optional character string naming a watershed column. If
#'   `NULL`, all rows are treated as one watershed.
#' @param na.rm Logical; remove missing values before calculation?
#' @param digits Non-negative integer giving the number of decimal places in
#'   calculated outputs.
#'
#' @return A list with `by_event`, containing event-level CVs and diagnostics,
#'   and `watershed_summary`, containing the unweighted mean and sample SD of
#'   eligible event CVs for each watershed.
#' @export
#'
#' @examples
#' x <- data.frame(
#'   event = c("a", "a", "a", "b", "b", "b"),
#'   concentration = c(1, 2, 3, 2, 3, 4)
#' )
#' spatial_cv(x, "concentration", "event")
spatial_cv <- function(
    data,
    concentration,
    event,
    watershed = NULL,
    na.rm = TRUE,
    digits = 2L
) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  validate_concentration_columns(concentration)
  validate_column_name(event, "event")
  if (!is.null(watershed)) {
    validate_column_name(watershed, "watershed")
  }
  validate_digits(digits)

  if (length(concentration) > 1L) {
    constituent_labels <- names(concentration)
    if (is.null(constituent_labels)) {
      constituent_labels <- concentration
    } else {
      unnamed <- is.na(constituent_labels) | !nzchar(constituent_labels)
      constituent_labels[unnamed] <- concentration[unnamed]
    }
    if (anyDuplicated(constituent_labels)) {
      stop("Constituent labels must be unique.", call. = FALSE)
    }

    results <- lapply(seq_along(concentration), function(i) {
      spatial_cv(
        data = data,
        concentration = unname(concentration[i]),
        event = event,
        watershed = watershed,
        na.rm = na.rm,
        digits = digits
      )
    })
    by_event <- do.call(rbind, Map(function(result, label) {
      data.frame(
        constituent = label,
        result$by_event,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }, results, constituent_labels))
    watershed_summary <- do.call(rbind, Map(function(result, label) {
      data.frame(
        constituent = label,
        result$watershed_summary,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }, results, constituent_labels))
    rownames(by_event) <- NULL
    rownames(watershed_summary) <- NULL

    return(structure(
      list(by_event = by_event, watershed_summary = watershed_summary),
      excluded_missing_group = attr(results[[1L]], "excluded_missing_group"),
      class = c("watershed_spatial_cv", "list")
    ))
  }

  required_columns <- c(concentration, event, watershed)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ", paste(missing_columns, collapse = ", "),
      ".", call. = FALSE
    )
  }
  if (!is.numeric(data[[concentration]])) {
    stop("The concentration column must be numeric.", call. = FALSE)
  }
  validate_na_rm(na.rm)

  valid_group <- !is.na(data[[event]])
  if (is.character(data[[event]]) || is.factor(data[[event]])) {
    valid_group <- valid_group & nzchar(trimws(as.character(data[[event]])))
  }
  if (!is.null(watershed)) {
    valid_group <- valid_group & !is.na(data[[watershed]])
    if (is.character(data[[watershed]]) || is.factor(data[[watershed]])) {
      valid_group <- valid_group &
        nzchar(trimws(as.character(data[[watershed]])))
    }
  }
  analysis_data <- data[valid_group, , drop = FALSE]
  if (nrow(analysis_data) == 0L) {
    stop("No rows have non-missing event and watershed identifiers.", call. = FALSE)
  }

  grouping_columns <- c(watershed, event)
  groups <- split(
    seq_len(nrow(analysis_data)),
    interaction(analysis_data[grouping_columns], drop = TRUE, lex.order = TRUE)
  )
  event_rows <- lapply(groups, function(index) {
    values <- analysis_data[[concentration]][index]
    n_missing <- sum(is.na(values))
    used <- if (na.rm) values[!is.na(values)] else values
    identifiers <- analysis_data[index[1L], grouping_columns, drop = FALSE]

    data.frame(
      identifiers,
      n_observations = length(values),
      n_used = if (na.rm) length(used) else sum(!is.na(used)),
      n_missing = n_missing,
      mean_concentration = round(
        if (length(used) == 0L) NA_real_ else mean(used), digits
      ),
      sd_concentration = round(
        if (length(used) < 2L) NA_real_ else stats::sd(used), digits
      ),
      .spatial_cv_unrounded = coefficient_of_variation(values, na.rm = na.rm),
      spatial_cv = round(
        coefficient_of_variation(values, na.rm = na.rm), digits
      ),
      check.names = FALSE
    )
  })
  by_event <- do.call(rbind, event_rows)
  rownames(by_event) <- NULL

  if (is.null(watershed)) {
    eligible <- is.finite(by_event$.spatial_cv_unrounded)
    watershed_summary <- data.frame(
      n_events = nrow(by_event),
      n_events_used = sum(eligible),
      mean_spatial_cv = if (any(eligible)) {
        round(mean(by_event$.spatial_cv_unrounded[eligible]), digits)
      } else {
        NA_real_
      },
      sd_spatial_cv = if (sum(eligible) >= 2L) {
        round(stats::sd(by_event$.spatial_cv_unrounded[eligible]), digits)
      } else {
        NA_real_
      }
    )
  } else {
    watershed_groups <- split(seq_len(nrow(by_event)), by_event[[watershed]])
    summaries <- lapply(watershed_groups, function(index) {
      eligible <- is.finite(by_event$.spatial_cv_unrounded[index])
      data.frame(
        by_event[index[1L], watershed, drop = FALSE],
        n_events = length(index),
        n_events_used = sum(eligible),
        mean_spatial_cv = if (any(eligible)) {
          round(mean(by_event$.spatial_cv_unrounded[index][eligible]), digits)
        } else {
          NA_real_
        },
        sd_spatial_cv = if (sum(eligible) >= 2L) {
          round(
            stats::sd(by_event$.spatial_cv_unrounded[index][eligible]),
            digits
          )
        } else {
          NA_real_
        },
        check.names = FALSE
      )
    })
    watershed_summary <- do.call(rbind, summaries)
    rownames(watershed_summary) <- NULL
  }
  by_event$.spatial_cv_unrounded <- NULL

  structure(
    list(by_event = by_event, watershed_summary = watershed_summary),
    excluded_missing_group = sum(!valid_group),
    class = c("watershed_spatial_cv", "list")
  )
}

#' Temporal coefficient of variation
#'
#' Calculates the coefficient of variation across sampling events separately
#' for each site using sample SD (denominator `n - 1`), then calculates the
#' unweighted mean of the eligible site CVs for each watershed.
#'
#' @param data A data frame containing concentration and site columns.
#' @param concentration Character string naming one numeric concentration
#'   column, or a named character vector mapping constituent labels to numeric
#'   concentration columns.
#' @param site Character string naming the sampling-site column.
#' @param watershed Optional character string naming a watershed column. If
#'   `NULL`, all rows are treated as one watershed.
#' @param na.rm Logical; remove missing values before calculation?
#' @param digits Non-negative integer giving the number of decimal places in
#'   calculated outputs.
#'
#' @return A list with `by_site`, containing site-level CVs and diagnostics,
#'   and `watershed_summary`, containing the unweighted mean and sample SD of
#'   eligible site CVs for each watershed.
#' @export
#'
#' @examples
#' x <- data.frame(
#'   site = c("a", "a", "a", "b", "b", "b"),
#'   concentration = c(1, 2, 3, 2, 3, 4)
#' )
#' temporal_cv(x, "concentration", "site")
temporal_cv <- function(
    data,
    concentration,
    site,
    watershed = NULL,
    na.rm = TRUE,
    digits = 2L
) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  validate_concentration_columns(concentration)
  validate_column_name(site, "site")
  if (!is.null(watershed)) {
    validate_column_name(watershed, "watershed")
  }
  validate_digits(digits)

  if (length(concentration) > 1L) {
    constituent_labels <- names(concentration)
    if (is.null(constituent_labels)) {
      constituent_labels <- concentration
    } else {
      unnamed <- is.na(constituent_labels) | !nzchar(constituent_labels)
      constituent_labels[unnamed] <- concentration[unnamed]
    }
    if (anyDuplicated(constituent_labels)) {
      stop("Constituent labels must be unique.", call. = FALSE)
    }

    results <- lapply(seq_along(concentration), function(i) {
      temporal_cv(
        data = data,
        concentration = unname(concentration[i]),
        site = site,
        watershed = watershed,
        na.rm = na.rm,
        digits = digits
      )
    })
    by_site <- do.call(rbind, Map(function(result, label) {
      data.frame(
        constituent = label,
        result$by_site,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }, results, constituent_labels))
    watershed_summary <- do.call(rbind, Map(function(result, label) {
      data.frame(
        constituent = label,
        result$watershed_summary,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }, results, constituent_labels))
    rownames(by_site) <- NULL
    rownames(watershed_summary) <- NULL

    return(structure(
      list(by_site = by_site, watershed_summary = watershed_summary),
      excluded_missing_group = attr(results[[1L]], "excluded_missing_group"),
      class = c("watershed_temporal_cv", "list")
    ))
  }

  required_columns <- c(concentration, site, watershed)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ", paste(missing_columns, collapse = ", "),
      ".", call. = FALSE
    )
  }
  if (!is.numeric(data[[concentration]])) {
    stop("The concentration column must be numeric.", call. = FALSE)
  }
  validate_na_rm(na.rm)

  valid_group <- !is.na(data[[site]])
  if (is.character(data[[site]]) || is.factor(data[[site]])) {
    valid_group <- valid_group & nzchar(trimws(as.character(data[[site]])))
  }
  if (!is.null(watershed)) {
    valid_group <- valid_group & !is.na(data[[watershed]])
    if (is.character(data[[watershed]]) || is.factor(data[[watershed]])) {
      valid_group <- valid_group &
        nzchar(trimws(as.character(data[[watershed]])))
    }
  }
  analysis_data <- data[valid_group, , drop = FALSE]
  if (nrow(analysis_data) == 0L) {
    stop("No rows have non-missing site and watershed identifiers.", call. = FALSE)
  }

  grouping_columns <- c(watershed, site)
  groups <- split(
    seq_len(nrow(analysis_data)),
    interaction(analysis_data[grouping_columns], drop = TRUE, lex.order = TRUE)
  )
  site_rows <- lapply(groups, function(index) {
    values <- analysis_data[[concentration]][index]
    n_missing <- sum(is.na(values))
    used <- if (na.rm) values[!is.na(values)] else values
    identifiers <- analysis_data[index[1L], grouping_columns, drop = FALSE]

    data.frame(
      identifiers,
      n_observations = length(values),
      n_used = if (na.rm) length(used) else sum(!is.na(used)),
      n_missing = n_missing,
      mean_concentration = round(
        if (length(used) == 0L) NA_real_ else mean(used), digits
      ),
      sd_concentration = round(
        if (length(used) < 2L) NA_real_ else stats::sd(used), digits
      ),
      .temporal_cv_unrounded = coefficient_of_variation(values, na.rm = na.rm),
      temporal_cv = round(
        coefficient_of_variation(values, na.rm = na.rm), digits
      ),
      check.names = FALSE
    )
  })
  by_site <- do.call(rbind, site_rows)
  rownames(by_site) <- NULL

  if (is.null(watershed)) {
    eligible <- is.finite(by_site$.temporal_cv_unrounded)
    watershed_summary <- data.frame(
      n_sites = nrow(by_site),
      n_sites_used = sum(eligible),
      mean_temporal_cv = if (any(eligible)) {
        round(mean(by_site$.temporal_cv_unrounded[eligible]), digits)
      } else {
        NA_real_
      },
      sd_temporal_cv = if (sum(eligible) >= 2L) {
        round(stats::sd(by_site$.temporal_cv_unrounded[eligible]), digits)
      } else {
        NA_real_
      }
    )
  } else {
    watershed_groups <- split(seq_len(nrow(by_site)), by_site[[watershed]])
    summaries <- lapply(watershed_groups, function(index) {
      eligible <- is.finite(by_site$.temporal_cv_unrounded[index])
      data.frame(
        by_site[index[1L], watershed, drop = FALSE],
        n_sites = length(index),
        n_sites_used = sum(eligible),
        mean_temporal_cv = if (any(eligible)) {
          round(mean(by_site$.temporal_cv_unrounded[index][eligible]), digits)
        } else {
          NA_real_
        },
        sd_temporal_cv = if (sum(eligible) >= 2L) {
          round(
            stats::sd(by_site$.temporal_cv_unrounded[index][eligible]),
            digits
          )
        } else {
          NA_real_
        },
        check.names = FALSE
      )
    })
    watershed_summary <- do.call(rbind, summaries)
    rownames(watershed_summary) <- NULL
  }
  by_site$.temporal_cv_unrounded <- NULL

  structure(
    list(by_site = by_site, watershed_summary = watershed_summary),
    excluded_missing_group = sum(!valid_group),
    class = c("watershed_temporal_cv", "list")
  )
}

coefficient_of_variation <- function(concentration, na.rm) {
  if (!is.numeric(concentration)) {
    stop("concentration must be numeric.", call. = FALSE)
  }
  validate_na_rm(na.rm)

  if (na.rm) {
    concentration <- concentration[!is.na(concentration)]
  }
  if (length(concentration) < 2L) {
    return(NA_real_)
  }

  concentration_mean <- mean(concentration)
  if (is.na(concentration_mean) || concentration_mean == 0) {
    return(NA_real_)
  }

  stats::sd(concentration) / concentration_mean
}

validate_column_name <- function(x, argument) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(argument, " must be one non-empty character string.", call. = FALSE)
  }
  invisible(x)
}

validate_concentration_columns <- function(x) {
  if (!is.character(x) || length(x) < 1L || anyNA(x) || any(!nzchar(x))) {
    stop(
      "concentration must contain one or more non-empty column names.",
      call. = FALSE
    )
  }
  if (anyDuplicated(unname(x))) {
    stop("Concentration columns must be unique.", call. = FALSE)
  }
  invisible(x)
}

validate_na_rm <- function(na.rm) {
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("na.rm must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(na.rm)
}

validate_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      digits < 0 || digits != as.integer(digits)) {
    stop("digits must be one non-negative integer.", call. = FALSE)
  }
  invisible(digits)
}
