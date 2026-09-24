#' Calculate subcatchment leverage from long-format data
#'
#' Expects a long-format data frame and calculates leverage separately within
#' every watershed-constituent-event
#' group and looks up that group's outlet observation. Inputs must contain one
#' row per watershed-constituent-event-site combination. Samples in an event
#' should represent comparable hydrologic conditions.
#'
#' `contributing_area_col` must contain total upstream drainage area, not
#' incremental area. `L` and `Lpct` use measured specific discharge; `LS` and
#' `LSpct` assume that specific discharge is equal between each site and its
#' catchment outlet. The function makes no source/sink classification.
#' The specific-discharge input must be calculated before calling this
#' function as \eqn{q = Q/A}, where \eqn{Q} is total volumetric discharge and
#' \eqn{A} is contributing area. The function does not convert total discharge
#' to specific discharge internally.
#'
#' Measured-specific-discharge leverage is
#' \deqn{L = (C_s-C_o)(A_s/A_o)(q_s/q_o).}
#' Here \eqn{q} is specific discharge (total discharge divided by contributing
#' area). Standardized leverage is
#' \deqn{LS = (C_s - C_o)(A_s / A_o).}
#' `Lpct` and `LSpct` divide the corresponding value by \eqn{C_o} and
#' multiply by 100. `L` and `LS` have concentration units; neither is solute
#' mass flux.
#'
#' @param data Long-format data frame.
#' @param concentration_col Name of the numeric concentration column.
#' @param site_col Name of the site column.
#' @param event_col Name of the sampling-event column. An event is one sampling
#'   campaign or synoptic sampling occasion.
#' @param watershed_col Name of the watershed column.
#' @param constituent_col Name of the constituent column.
#' @param contributing_area_col Name of the numeric total upstream drainage-area
#'   column. Site and outlet area units must match.
#' @param outlet_site One outlet identifier used in all watersheds, or a vector
#'   of outlet identifiers named by watershed.
#' @param variant One of `"L"`, `"Lpct"`, `"LS"`, or `"LSpct"`.
#' @param specific_discharge_col Name of the numeric specific-discharge column.
#'   Values must already be calculated as total volumetric discharge divided by
#'   contributing area (\eqn{q=Q/A}); do not supply unstandardized total
#'   discharge. Required for `L` and `Lpct`; unused for `LS` and `LSpct`.
#'
#' @return `data` with outlet values, `leverage_variant`, and `leverage`
#'   appended. `outlet_specific_discharge` is `NA` for standardized variants.
#' @export
#'
#' @examples
#' d <- data.frame(
#'   watershed = "w1", constituent = "DOC", event = "e1",
#'   site = c("a", "b", "outlet"), concentration = c(2, 3, 2.5),
#'   contributing_area = c(10, 20, 100), specific_discharge = c(0.8, 1.2, 1)
#' )
#' subcatchment_leverage(
#'   d, "concentration", "site", "event", "watershed", "constituent",
#'   "contributing_area", outlet_site = "outlet", variant = "LSpct"
#' )
subcatchment_leverage <- function(
    data,
    concentration_col,
    site_col,
    event_col,
    watershed_col,
    constituent_col,
    contributing_area_col,
    outlet_site,
    variant = c("LSpct", "LS", "L", "Lpct"),
    specific_discharge_col = NULL
) {
  variant <- match.arg(variant)
  uses_specific_discharge <- variant %in% c("L", "Lpct")
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (uses_specific_discharge &&
      (is.null(specific_discharge_col) || length(specific_discharge_col) != 1L ||
       is.na(specific_discharge_col) || specific_discharge_col == "")) {
    stop("specific_discharge_col is required for variant '", variant, "'.",
         call. = FALSE)
  }

  column_args <- c(
    concentration_col = concentration_col,
    site_col = site_col,
    event_col = event_col,
    watershed_col = watershed_col,
    constituent_col = constituent_col,
    contributing_area_col = contributing_area_col
  )
  if (uses_specific_discharge) {
    column_args <- c(
      column_args,
      specific_discharge_col = specific_discharge_col
    )
  }
  if (any(lengths(column_args) != 1L) || any(is.na(column_args)) ||
      any(column_args == "")) {
    stop("Each required column argument must be one non-missing column name.",
         call. = FALSE)
  }
  missing_cols <- setdiff(unname(column_args), names(data))
  if (length(missing_cols)) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  concentration <- data[[concentration_col]]
  site <- data[[site_col]]
  event <- data[[event_col]]
  watershed <- data[[watershed_col]]
  constituent <- data[[constituent_col]]
  area <- data[[contributing_area_col]]
  specific_discharge <- if (uses_specific_discharge) {
    data[[specific_discharge_col]]
  } else {
    NULL
  }

  numeric_values <- list(concentration = concentration, contributing_area = area)
  if (uses_specific_discharge) {
    numeric_values$specific_discharge <- specific_discharge
  }
  if (!all(vapply(numeric_values, is.numeric, logical(1)))) {
    stop("Concentration, area, and required specific-discharge columns must be numeric.",
         call. = FALSE)
  }
  if (any(area <= 0, na.rm = TRUE)) {
    stop("Contributing area must be positive.", call. = FALSE)
  }
  if (uses_specific_discharge && any(specific_discharge < 0, na.rm = TRUE)) {
    stop("Specific discharge cannot be negative.", call. = FALSE)
  }

  identifiers <- list(
    watershed = watershed, constituent = constituent, event = event, site = site
  )
  missing_identifier <- Reduce(`|`, lapply(identifiers, function(x) {
    is.na(x) | !nzchar(trimws(as.character(x)))
  }))
  if (any(missing_identifier)) {
    stop("Watershed, constituent, event, and site cannot be missing or blank.",
         call. = FALSE)
  }
  observation_key <- do.call(
    interaction, c(identifiers, list(drop = TRUE, lex.order = TRUE))
  )
  if (anyDuplicated(observation_key)) {
    duplicate_rows <- which(duplicated(observation_key) |
                              duplicated(observation_key, fromLast = TRUE))
    stop("Each watershed-constituent-event-site combination must be unique; ",
         "duplicate rows: ", paste(duplicate_rows, collapse = ", "), ".",
         call. = FALSE)
  }

  watersheds <- unique(as.character(watershed))
  if (length(outlet_site) == 1L && !is.na(outlet_site)) {
    outlet_by_watershed <- setNames(rep(outlet_site, length(watersheds)), watersheds)
  } else {
    if (is.null(names(outlet_site)) || any(names(outlet_site) == "") ||
        anyDuplicated(names(outlet_site))) {
      stop("Multiple outlet_site values must be uniquely named by watershed.",
           call. = FALSE)
    }
    missing_outlets <- setdiff(watersheds, names(outlet_site))
    if (length(missing_outlets)) {
      stop("No outlet supplied for watersheds: ",
           paste(missing_outlets, collapse = ", "), ".", call. = FALSE)
    }
    outlet_by_watershed <- outlet_site[watersheds]
  }

  group_key <- interaction(watershed, constituent, event,
                           drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), group_key)
  result <- data
  result$outlet_site <- unname(outlet_by_watershed[as.character(watershed)])
  result$outlet_concentration <- NA_real_
  result$outlet_area <- NA_real_
  result$outlet_specific_discharge <- NA_real_

  for (rows in groups) {
    target <- result$outlet_site[rows[1L]]
    outlet_rows <- rows[as.character(site[rows]) == as.character(target)]
    label <- paste0("watershed '", watershed[rows[1L]], "', constituent '",
                    constituent[rows[1L]], "', event '", event[rows[1L]], "'")
    if (length(outlet_rows) != 1L) {
      stop("Group ", label, " must contain exactly one outlet row; found ",
           length(outlet_rows), ".", call. = FALSE)
    }
    outlet_row <- outlet_rows[[1L]]
    incomplete <- is.na(concentration[outlet_row]) || is.na(area[outlet_row]) ||
      (uses_specific_discharge && is.na(specific_discharge[outlet_row]))
    if (incomplete) {
      stop("The outlet row for ", label, " has a missing required value.",
           call. = FALSE)
    }
    if (uses_specific_discharge && specific_discharge[outlet_row] <= 0) {
      stop("Outlet specific discharge must be positive for ", label, ".",
           call. = FALSE)
    }
    if (variant %in% c("Lpct", "LSpct") && concentration[outlet_row] == 0) {
      stop("Outlet concentration cannot be zero for percent leverage in ",
           label, ".", call. = FALSE)
    }
    result$outlet_concentration[rows] <- concentration[outlet_row]
    result$outlet_area[rows] <- area[outlet_row]
    if (uses_specific_discharge) {
      result$outlet_specific_discharge[rows] <- specific_discharge[outlet_row]
    }
  }

  weight <- if (uses_specific_discharge) {
    (area / result$outlet_area) *
      (specific_discharge / result$outlet_specific_discharge)
  } else {
    area / result$outlet_area
  }
  result$leverage_variant <- rep(variant, nrow(data))
  result$leverage <- (concentration - result$outlet_concentration) * weight
  if (variant %in% c("Lpct", "LSpct")) {
    result$leverage <- 100 * result$leverage / result$outlet_concentration
  }
  result
}
