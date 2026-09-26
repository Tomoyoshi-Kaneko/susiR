# ============================================================================
# raw_import.R -- convert an arbitrary plate-reader export (unknown sheet
# name, unknown header row, unknown well-ID spelling, unknown time column
# position/units, extra metadata columns before/after the data) into the
# standard two-sheet shape (R = wide time-course, name = well mapping
# template) that read_plate_data() expects. This is the layer that lets
# susiR accept files from readers other than the one used so far, without
# needing every downstream function to know about reader-specific layouts.
#
# Scope, deliberately: this handles the *OD time-course table* only. The
# well -> condition mapping ("name" sheet) is never something plate-reader
# software produces -- it always has to come from the experimenter -- so
# this writes an empty template (well IDs filled in, sample/bio_rep blank)
# for the user to complete, rather than trying to guess it.
#
# Reads via readxl when available, falling back to openxlsx: real-world
# exports (large sheets, embedded control characters in headers, sparse
# columns) were found to trip up openxlsx's C parser on at least one real
# file in testing, while readxl handled it without issue.
# ============================================================================

#' @keywords internal
.excel_col_letter <- function(idx) {
  letters_out <- character(length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]; s <- ""
    while (i > 0) {
      rem <- (i - 1) %% 26
      s <- paste0(LETTERS[rem + 1], s)
      i <- (i - 1) %/% 26
    }
    letters_out[k] <- s
  }
  letters_out
}

#' @keywords internal
.sheet_names_any <- function(file_path) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    out <- tryCatch(readxl::excel_sheets(file_path), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  openxlsx::getSheetNames(file_path)
}

#' Read a sheet as a raw character/mixed data frame, no header assumptions
#'
#' Tries `readxl` first (found to be more robust on messy real-world plate
#' reader exports during development), falls back to `openxlsx`.
#' @keywords internal
.read_sheet_raw <- function(file_path, sheet) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    out <- tryCatch(
      as.data.frame(readxl::read_excel(file_path, sheet = sheet, col_names = FALSE,
                                        col_types = "text", .name_repair = "minimal"),
                     stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }
  tryCatch(
    openxlsx::read.xlsx(file_path, sheet = sheet, colNames = FALSE, skipEmptyRows = FALSE),
    error = function(e) NULL
  )
}

#' Normalize a well-ID-like string to the canonical "A1" style
#'
#' Accepts plain well IDs ("A1", "A01", "A:1", "a 1") as well as a well ID
#' embedded inside a longer label, e.g. a parenthesized form such as
#' `"Sample0001 (B02)"` or `"\u30b5\u30f3\u30d7\u30eb0001 (B02)"` (as produced by at
#' least one real Thermo Fisher SkanIt export) -- the parenthesized form is
#' tried first since the surrounding text (which may itself contain
#' digits, as here) would otherwise confuse a whole-string match. Returns
#' `NA` for strings that don't look like a well ID at all (row letter
#' A-P, column 1-24).
#'
#' @param x Character vector.
#' @return Character vector, same length, normalized well IDs (or `NA`).
#' @export
susi_normalize_well_id <- function(x) {
  one <- function(s) {
    if (is.na(s)) return(NA_character_)
    paren <- regmatches(s, regexpr("\\(([^()]*)\\)", s))
    if (length(paren) == 1 && nzchar(paren)) {
      inner <- toupper(gsub("[^A-Pa-p0-9]", "", substring(paren, 2, nchar(paren) - 1)))
      mm <- regmatches(inner, regexec("^([A-P])0*([0-9]{1,2})$", inner))[[1]]
      if (length(mm) == 3) return(paste0(mm[2], mm[3]))
    }
    s2 <- toupper(trimws(s))
    s2 <- gsub("[^A-P0-9]", "", s2)
    mm2 <- regmatches(s2, regexec("^([A-P])0*([0-9]{1,2})$", s2))[[1]]
    if (length(mm2) == 3) return(paste0(mm2[2], mm2[3]))
    NA_character_
  }
  vapply(as.character(x), one, character(1), USE.NAMES = FALSE)
}

#' Find a plate-reader data table inside a workbook of unknown shape
#'
#' Scans every sheet (or just `sheet` if given) for a row that looks like a
#' time-course header: contains the word "time" in some cell, and has at
#' least `min_wells` cells that look like well IDs ([susi_normalize_well_id()]
#' returns non-`NA`). This is intentionally conservative -- a low-confidence
#' or ambiguous match is reported as not-found rather than guessed at,
#' since silently parsing the wrong table would be worse than stopping.
#'
#' @param file_path Path to the `.xlsx` workbook.
#' @param sheet Optional sheet name to restrict the search to.
#' @param min_wells Minimum number of well-ID-like header cells required to
#'   accept a row as the header row. Default 6.
#' @return A list (`sheet`, `header_row`, `n_wells`) on success, or `NULL`
#'   if nothing confident enough was found.
#' @export
susi_find_data_table <- function(file_path, sheet = NULL, min_wells = 6) {
  candidates <- if (!is.null(sheet)) sheet else .sheet_names_any(file_path)
  for (sn in candidates) {
    raw <- .read_sheet_raw(file_path, sn)
    if (is.null(raw) || nrow(raw) == 0) next
    n_scan <- min(nrow(raw), 150)
    for (r in seq_len(n_scan)) {
      row_vals <- as.character(unlist(raw[r, ]))
      has_time <- any(grepl("time|\u6642\u9593", row_vals, ignore.case = TRUE))
      n_wells <- sum(!is.na(susi_normalize_well_id(row_vals)))
      if (has_time && n_wells >= min_wells) {
        return(list(sheet = sn, header_row = r, n_wells = n_wells))
      }
    }
  }
  NULL
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Import a plate-reader export of unknown format
#'
#' Locates the data table ([susi_find_data_table()]), identifies which
#' columns are wells (any spelling [susi_normalize_well_id()] recognizes),
#' which column is time (header text containing "time"; unit guessed from
#' "sec"/"min"/"h" in that header text if present, otherwise from the
#' workbook's actual Excel cell formatting via [susi_resolve_time_seconds()]),
#' and drops everything else (temperature columns, read-type/read-number
#' columns, trailing metadata) -- then writes a workbook in the standard
#' shape: an `"R"` sheet (`Time` in seconds + one column per well) plus a
#' `"name"` sheet template (well IDs filled in; `sample`/`bio_rep` columns
#' present but blank, for the user to complete). The result opens directly
#' with [read_plate_data()]/[run_susi()] -- no further changes to this
#' package are needed for a new reader brand, only this one-time conversion
#' step per raw file.
#'
#' Always inspect the printed summary (sheet/row found, well count, time
#' range) before trusting the output -- for an unfamiliar or unusually
#' laid-out export, auto-detection can fail loudly (a clear error) but
#' could in principle also mis-detect quietly, and the summary is the
#' check against that.
#'
#' @param file_path Path to the raw `.xlsx` export.
#' @param output_path Path to write the converted workbook to.
#' @param sheet Optional: restrict the search to one sheet (use when the
#'   file has several sheets and auto-detection finds the wrong one, or to
#'   skip the scan when you already know which sheet holds the data).
#' @param time_unit Optional override: `"seconds"`, `"minutes"`, or
#'   `"hours"`. If not given, guessed from the time column's header text,
#'   falling back to automatic format-based detection.
#' @param overwrite Overwrite `output_path` if it already exists.
#' @return (Invisibly) a list: `wells` (normalized well IDs found),
#'   `time_range_hours`, `source_sheet`, `header_row`, `output_path`. Also
#'   prints a human-readable summary.
#' @export
susi_import_raw <- function(file_path, output_path, sheet = NULL, time_unit = NULL, overwrite = TRUE) {
  found <- susi_find_data_table(file_path, sheet = sheet)
  if (is.null(found)) {
    stop("Could not find a plate-reader time-course table in this file",
         if (!is.null(sheet)) paste0(" (sheet '", sheet, "')") else " (scanned all sheets)",
         ". Looked for a row containing a 'time' label and at least 6 well-ID-like ",
         "column headers (e.g. A1, A:1, A01). If this file uses an unusual layout, ",
         "either pass `sheet` to point at the right tab, or convert it manually to ",
         "match the T1_rawdata.xlsx example format.", call. = FALSE)
  }

  raw <- .read_sheet_raw(file_path, found$sheet)
  header <- as.character(unlist(raw[found$header_row, ]))
  header_clean <- gsub("[\r\n]+", " ", header)
  body <- raw[seq(found$header_row + 1, nrow(raw)), , drop = FALSE]
  names(body) <- header

  well_ids <- susi_normalize_well_id(header)
  is_well_col <- !is.na(well_ids)
  time_col_idx <- which(grepl("time|\u6642\u9593", header_clean, ignore.case = TRUE))[1]
  if (is.na(time_col_idx)) {
    stop("Found ", sum(is_well_col), " well-like columns in sheet '", found$sheet,
         "' (row ", found$header_row, ") but no column header containing 'time' (or \u6642\u9593).", call. = FALSE)
  }
  time_header_text <- trimws(header_clean[time_col_idx])
  time_raw <- suppressWarnings(as.numeric(body[[time_col_idx]]))

  unit <- time_unit
  method <- "explicit_override"
  if (is.null(unit)) {
    if (grepl("sec|\u79d2", time_header_text, ignore.case = TRUE)) { unit <- "seconds"; method <- "header_text" }
    else if (grepl("min|\u5206", time_header_text, ignore.case = TRUE)) { unit <- "minutes"; method <- "header_text" }
    else if (grepl("^h(our)?s?$|\\(h\\)|\\bhr|\u6642\u9593\\s*\\[?h", time_header_text, ignore.case = TRUE)) { unit <- "hours"; method <- "header_text" }
  }
  if (is.null(unit)) {
    col_letter <- .excel_col_letter(time_col_idx)
    tres <- susi_resolve_time_seconds(time_raw, file_path = file_path, sheet = found$sheet, col_letter = col_letter)
    time_sec <- tres$time_sec
    method <- tres$method
  } else {
    time_sec <- switch(unit, seconds = time_raw, minutes = time_raw * 60, hours = time_raw * 3600)
  }

  well_data <- body[, is_well_col, drop = FALSE]
  names(well_data) <- well_ids[is_well_col]
  well_data <- as.data.frame(lapply(well_data, function(x) suppressWarnings(as.numeric(x))))

  if (anyDuplicated(names(well_data))) {
    dups <- unique(names(well_data)[duplicated(names(well_data))])
    warning("Duplicate well IDs after normalization (kept first occurrence): ", paste(dups, collapse = ", "), call. = FALSE)
    well_data <- well_data[, !duplicated(names(well_data)), drop = FALSE]
  }

  ok_rows <- !is.na(time_sec)
  out_od <- cbind(data.frame(Time = time_sec[ok_rows]), well_data[ok_rows, , drop = FALSE])

  name_template <- data.frame(well = names(well_data), sample = "", bio_rep = "", stringsAsFactors = FALSE)

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "R")
  openxlsx::writeData(wb, "R", out_od)
  openxlsx::addWorksheet(wb, "name")
  openxlsx::writeData(wb, "name", name_template)
  openxlsx::saveWorkbook(wb, output_path, overwrite = overwrite)

  time_range_h <- range(time_sec[ok_rows]) / 3600
  message(sprintf(
    "Imported from sheet '%s' (header row %d): %d wells, time column '%s' (unit: %s, detected via %s), %.2f to %.2f hours over %d points.\nWrote %s -- fill in the 'sample' (and optionally 'bio_rep') column in its 'name' sheet, then use it with read_plate_data()/run_susi() as normal.",
    found$sheet, found$header_row, ncol(well_data), time_header_text, unit %||% "auto", method,
    time_range_h[1], time_range_h[2], sum(ok_rows), output_path
  ))

  invisible(list(wells = names(well_data), time_range_hours = time_range_h,
                  source_sheet = found$sheet, header_row = found$header_row, output_path = output_path))
}
