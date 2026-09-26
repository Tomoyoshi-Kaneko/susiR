# ============================================================================
# data_input.R -- reading + validating plate data, generalized replicate
#                 resolution (explicit bio_rep column, with a positional
#                 fallback that mirrors the bundled example dataset's layout)
# ============================================================================

#' Read a plate-reader OD time course together with its sample-mapping table
#'
#' The workbook must contain two sheets: a wide-format OD time course
#' (`od_sheet`; first column = time, every other column = one well) and a
#' sample-mapping table (`map_sheet`) with one row per well. The mapping
#' table needs, at minimum, a well-ID column (`well_col`) and a
#' condition-label column (`condition_col`); condition labels are treated as
#' arbitrary strings (MOI values, phage names, anything) and are never
#' assumed to belong to a fixed set.
#'
#' Biological-replicate membership is resolved by [susi_resolve_bio_rep()]:
#' if `bio_rep_col` is present and filled in, it is used directly (allowing
#' any number of replicates, laid out however the user likes); otherwise
#' replicate membership is inferred using the same convention as the bundled
#' example dataset (consecutive blocks of `tech_reps_per_bio_rep` wells, in
#' plate order, per condition) -- see [susi_resolve_bio_rep()] for details.
#'
#' @param file_path Path to the `.xlsx` workbook.
#' @param od_sheet Name of the OD time-course sheet. Default `"R"`.
#' @param map_sheet Name of the sample-mapping sheet. Default `"name"`.
#' @param well_col Column in `map_sheet` holding well IDs. Default `"well"`.
#' @param condition_col Column in `map_sheet` holding condition/sample
#'   labels. Default `"sample"`.
#' @param bio_rep_col Column in `map_sheet` holding an explicit
#'   biological-replicate identifier, if present. Default `"bio_rep"`.
#' @param tech_reps_per_bio_rep Block size used for the positional fallback
#'   when `bio_rep_col` is absent/empty. Default `12`, matching the bundled
#'   example dataset (3 biological x 12 technical replicates per condition).
#' @param time_unit `"auto"` (default) inspects the workbook's actual Excel
#'   cell formatting to tell whether the Time column is a plain number
#'   (assumed to be elapsed seconds) or a genuine Excel time/duration cell
#'   (whose raw value is a fraction of a day), converting correctly either
#'   way -- see [susi_resolve_time_seconds()]. Set explicitly to
#'   `"seconds"` or `"hours"` to skip detection and force an interpretation.
#' @return A list with elements:
#'   \item{od_long}{Tidy long data frame: `well`, `time_h`, `od`, `condition`, `bio_rep`.}
#'   \item{mapping}{The resolved mapping table (one row per well, with a
#'     resolved `bio_rep` column).}
#'   \item{time_h}{Sorted unique observation time points, in hours.}
#' @export
read_plate_data <- function(file_path,
                             od_sheet = "R",
                             map_sheet = "name",
                             well_col = "well",
                             condition_col = "sample",
                             bio_rep_col = "bio_rep",
                             tech_reps_per_bio_rep = 12,
                             time_unit = c("auto", "seconds", "hours")) {
  time_unit <- match.arg(time_unit)
  if (!file.exists(file_path)) stop("File not found: ", file_path, call. = FALSE)

  od_wide <- openxlsx::read.xlsx(file_path, sheet = od_sheet)
  mapping <- openxlsx::read.xlsx(file_path, sheet = map_sheet)

  if (!"Time" %in% names(od_wide)) {
    stop("OD sheet '", od_sheet, "' must have a 'Time' column as its first column.", call. = FALSE)
  }
  if (!well_col %in% names(mapping)) {
    stop("Mapping sheet '", map_sheet, "' has no column '", well_col,
         "'. Set `well_col` to the correct column name.", call. = FALSE)
  }
  if (!condition_col %in% names(mapping)) {
    stop("Mapping sheet '", map_sheet, "' has no column '", condition_col,
         "'. Set `condition_col` to the correct column name.", call. = FALSE)
  }

  wells_in_od  <- setdiff(names(od_wide), "Time")
  wells_in_map <- as.character(mapping[[well_col]])
  missing_from_od  <- setdiff(wells_in_map, wells_in_od)
  missing_from_map <- setdiff(wells_in_od, wells_in_map)
  if (length(missing_from_od) > 0) {
    warning(length(missing_from_od), " well(s) listed in the mapping sheet have no OD data and will be dropped: ",
            paste(utils::head(missing_from_od, 10), collapse = ", "),
            if (length(missing_from_od) > 10) ", ..." else "", call. = FALSE)
  }
  if (length(missing_from_map) > 0) {
    message(length(missing_from_map), " well(s) have OD data but no entry in the mapping sheet and will be ignored: ",
            paste(utils::head(missing_from_map, 10), collapse = ", "),
            if (length(missing_from_map) > 10) ", ..." else "")
  }

  mapping <- mapping[wells_in_map %in% wells_in_od, , drop = FALSE]
  mapping[[well_col]] <- as.character(mapping[[well_col]])
  mapping[[condition_col]] <- as.character(mapping[[condition_col]])

  mapping <- susi_resolve_bio_rep(mapping,
                                   well_col = well_col,
                                   condition_col = condition_col,
                                   bio_rep_col = bio_rep_col,
                                   tech_reps_per_bio_rep = tech_reps_per_bio_rep)

  time_raw <- suppressWarnings(as.numeric(od_wide[["Time"]]))
  ## Some real-world exports append trailing summary rows below the actual
  ## kinetic data (e.g. software-computed "Max V", "R-Squared", "Lagtime"
  ## statistics, sometimes with plate-row letters as row labels). Any row
  ## whose Time value isn't a number can't be a real time point, so such
  ## rows are dropped here rather than treated as a parse failure -- this
  ## also naturally discards fully blank trailing rows.
  junk_rows <- is.na(time_raw)
  if (any(junk_rows)) {
    message(sum(junk_rows), " row(s) with a non-numeric Time value were dropped ",
            "(commonly trailing summary rows some plate-reader software appends below the data): ",
            paste(utils::head(as.character(od_wide[["Time"]][junk_rows]), 5), collapse = ", "),
            if (sum(junk_rows) > 5) ", ..." else "")
    od_wide <- od_wide[!junk_rows, , drop = FALSE]
    time_raw <- time_raw[!junk_rows]
  }

  ## Well columns occasionally contain a non-numeric marker for an invalid
  ## read (e.g. sensor saturation/overflow, often shown as "?????" or
  ## similar) instead of a number. These become NA (a missing OD reading
  ## for that well/time point) rather than corrupting the column's type or
  ## aborting the whole import.
  well_cols_all <- setdiff(names(od_wide), "Time")
  na_introduced <- 0L
  for (wc in well_cols_all) {
    orig <- od_wide[[wc]]
    if (!is.numeric(orig)) {
      converted <- suppressWarnings(as.numeric(orig))
      na_introduced <- na_introduced + sum(is.na(converted) & !is.na(orig) & trimws(as.character(orig)) != "")
      od_wide[[wc]] <- converted
    }
  }
  if (na_introduced > 0) {
    message(na_introduced, " well reading(s) were non-numeric (e.g. an instrument error/overflow marker) ",
            "and were treated as missing (NA) rather than aborting.")
  }

  if (time_unit == "auto") {
    time_res <- susi_resolve_time_seconds(time_raw, file_path = file_path, sheet = od_sheet, col_letter = "A")
    time_sec <- time_res$time_sec
    message(sprintf("Time column interpreted via %s: %.2f to %.2f hours over %d points.",
                     time_res$method, min(time_sec) / 3600, max(time_sec) / 3600, length(unique(time_sec))))
  } else {
    time_sec <- if (time_unit == "seconds") time_raw else time_raw * 3600
  }
  time_h <- time_sec / 3600
  od_wide[["Time"]] <- time_raw  # normalize to numeric so pivoting/lookup stay consistent
  # Lookup from each raw Time value to its resolved seconds, reused below so
  # the (potentially expensive) XML introspection only runs once per file.
  time_lookup <- stats::setNames(time_sec, as.character(time_raw))

  keep_wells <- mapping[[well_col]]
  od_long <- tidyr::pivot_longer(
    od_wide[, c("Time", intersect(keep_wells, names(od_wide))), drop = FALSE],
    cols = -"Time", names_to = "well", values_to = "od"
  )
  od_long$time_h <- unname(time_lookup[as.character(od_long$Time)]) / 3600
  od_long$Time <- NULL

  od_long <- dplyr::left_join(
    od_long,
    mapping[, c(well_col, condition_col, ".bio_rep")],
    by = stats::setNames(well_col, "well")
  )
  names(od_long)[names(od_long) == condition_col] <- "condition"
  names(od_long)[names(od_long) == ".bio_rep"] <- "bio_rep"

  list(
    od_long = od_long[, c("well", "time_h", "od", "condition", "bio_rep")],
    mapping = mapping,
    time_h = sort(unique(time_h))
  )
}

#' Resolve biological-replicate membership for each well
#'
#' If `bio_rep_col` exists in `mapping` and is filled in for every well
#' (within at least one condition), those values are used as-is -- this is
#' how a user declares an arbitrary number of biological replicates (1, 2,
#' 4, ...) laid out however they like on the plate. Any condition for which
#' the column is entirely missing/blank instead falls back to a positional
#' rule that mirrors the bundled example dataset: wells belonging to that
#' condition are ordered in plate order (row letter, then column number) and
#' split into consecutive blocks of `tech_reps_per_bio_rep` wells, numbered
#' 1, 2, 3, .... A condition with 36 wells and the default block size of 12
#' therefore yields 3 biological replicates (as in the example dataset); 12
#' wells yields 1; 24 wells yields 2; 48 wells yields 4; and so on. Mixing
#' explicit and positional resolution *within a single condition* is not
#' allowed (it is ambiguous) and raises an error asking the user to either
#' fill in `bio_rep` for every well of that condition or remove it entirely.
#'
#' @inheritParams read_plate_data
#' @param mapping Sample-mapping data frame (as read from the workbook).
#' @return `mapping` with an added `.bio_rep` integer column.
#' @export
susi_resolve_bio_rep <- function(mapping, well_col = "well", condition_col = "sample",
                                  bio_rep_col = "bio_rep", tech_reps_per_bio_rep = 12) {
  has_col <- bio_rep_col %in% names(mapping)
  explicit_vals <- if (has_col) mapping[[bio_rep_col]] else rep(NA, nrow(mapping))

  mapping$.bio_rep <- NA_integer_
  conditions <- unique(mapping[[condition_col]])

  fallback_used <- character(0)
  explicit_used <- character(0)

  for (cond in conditions) {
    idx <- which(mapping[[condition_col]] == cond)
    vals <- explicit_vals[idx]
    any_given <- any(!is.na(vals) & vals != "")
    all_given <- all(!is.na(vals) & vals != "")

    if (any_given && !all_given) {
      missing_wells <- mapping[[well_col]][idx][is.na(vals) | vals == ""]
      stop("Condition '", cond, "' has '", bio_rep_col, "' filled in for some wells but not others ",
           "(missing for: ", paste(utils::head(missing_wells, 10), collapse = ", "), "). ",
           "Either fill in '", bio_rep_col, "' for every well of this condition, or leave it blank for all of them ",
           "so the positional fallback (blocks of ", tech_reps_per_bio_rep, ") is used.", call. = FALSE)
    }

    if (all_given) {
      mapping$.bio_rep[idx] <- as.integer(vals)
      explicit_used <- c(explicit_used, cond)
    } else {
      ord <- idx[susi_well_sort_key(mapping[[well_col]][idx])]
      n <- length(ord)
      block <- pmin(ceiling(seq_len(n) / tech_reps_per_bio_rep), ceiling(n / tech_reps_per_bio_rep))
      mapping$.bio_rep[ord] <- ceiling(seq_len(n) / tech_reps_per_bio_rep)
      fallback_used <- c(fallback_used, cond)
    }
  }

  if (length(fallback_used) > 0) {
    n_reps <- sapply(fallback_used, function(cond) max(mapping$.bio_rep[mapping[[condition_col]] == cond]))
    message("Biological replicates inferred positionally (blocks of ", tech_reps_per_bio_rep,
            " wells, plate order) for: ",
            paste(sprintf("%s (n=%d)", fallback_used, n_reps), collapse = ", "))
  }
  if (length(explicit_used) > 0) {
    message("Biological replicates taken directly from column '", bio_rep_col, "' for: ",
            paste(explicit_used, collapse = ", "))
  }

  mapping
}
