# ============================================================================
# time_format.R -- robust handling of the "Time" column's units.
#
# Motivation: openxlsx::read.xlsx() returns whatever raw numeric value is
# stored in a cell, with NO unit conversion. Two real files in this
# project's own "already-correct" R-sheet format store semantically
# identical timestamps completely differently:
#   - T1_rawdata.xlsx:      Time cells are plain numbers, "General" format
#                            -> raw value IS elapsed seconds (e.g. 900)
#   - another lab's export: Time cells are genuine Excel duration cells,
#                            format code "[h]:mm" -> raw value is a
#                            fraction of a day (900 sec -> 0.010417)
# Treating the second case as if it were the first silently produces time
# values wrong by a factor of 86400 -- not an error, just silently corrupt
# results. This file detects which case applies from the workbook's actual
# XML formatting metadata (the authoritative source), with a magnitude
# heuristic as a fallback for when that can't be read.
# ============================================================================

.excel_builtin_time_formats <- c(
  `18` = "h:mm AM/PM", `19` = "h:mm:ss AM/PM", `20` = "h:mm", `21` = "h:mm:ss",
  `22` = "m/d/yy h:mm", `45` = "mm:ss", `46` = "[h]:mm:ss", `47` = "mmss.0"
)

#' Does an Excel number-format code represent a time/duration (not a date)?
#' @keywords internal
.is_time_format_code <- function(code) {
  if (is.na(code) || !nzchar(code)) return(FALSE)
  # Duration/time formats always include a literal ':' between h/m/s tokens
  # (distinguishing them from pure date formats like "mm-dd-yy") and use
  # h/s tokens (or the elapsed-time bracket form "[h]") somewhere.
  grepl(":", code, fixed = TRUE) && grepl("[hHsS]", code)
}

#' Read a sheet's per-cell number-format codes for one column, from the
#' workbook's raw XML (styles.xml + the sheet XML), without needing tidyxl.
#'
#' @param file_path Path to the .xlsx file.
#' @param sheet Sheet name.
#' @param col_letter Column letter (e.g. "A") whose format to check.
#' @param n_check How many data rows (after the header) to sample.
#' @return Character vector of format codes found (may be empty if the
#'   workbook can't be introspected this way -- callers should fall back).
#' @keywords internal
susi_excel_cell_formats <- function(file_path, sheet, col_letter = "A", n_check = 5) {
  result <- tryCatch({
    tmp <- tempfile(); dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    utils::unzip(file_path, exdir = tmp)

    wb_xml <- xml2::read_xml(file.path(tmp, "xl", "workbook.xml"))
    ns <- xml2::xml_ns(wb_xml)
    sheets <- xml2::xml_find_all(wb_xml, ".//d1:sheets/d1:sheet", ns)
    sheet_names <- xml2::xml_attr(sheets, "name")
    sheet_idx <- match(sheet, sheet_names)
    if (is.na(sheet_idx)) return(character(0))

    sheet_files <- list.files(file.path(tmp, "xl", "worksheets"), pattern = "^sheet[0-9]+\\.xml$", full.names = TRUE)
    sheet_files <- sheet_files[order(as.integer(gsub("\\D", "", basename(sheet_files))))]
    if (sheet_idx > length(sheet_files)) return(character(0))

    styles_xml <- xml2::read_xml(file.path(tmp, "xl", "styles.xml"))
    fmt_map <- as.list(.excel_builtin_time_formats)
    custom <- xml2::xml_find_all(styles_xml, ".//d1:numFmts/d1:numFmt", ns)
    for (nf in custom) fmt_map[[xml2::xml_attr(nf, "numFmtId")]] <- xml2::xml_attr(nf, "formatCode")
    xfs <- xml2::xml_find_all(styles_xml, ".//d1:cellXfs/d1:xf", ns)
    xf_numfmt <- xml2::xml_attr(xfs, "numFmtId")

    sdoc <- xml2::read_xml(sheet_files[sheet_idx])
    rows <- xml2::xml_find_all(sdoc, ".//d1:row", ns)
    codes <- character(0)
    checked <- 0
    for (row in rows) {
      cells <- xml2::xml_find_all(row, "d1:c", ns)
      refs <- xml2::xml_attr(cells, "r")
      hit <- which(grepl(paste0("^", col_letter, "[0-9]+$"), refs))
      for (h in hit) {
        s <- xml2::xml_attr(cells[[h]], "s")
        if (!is.na(s)) {
          idx <- as.integer(s) + 1
          if (idx <= length(xf_numfmt)) {
            code <- fmt_map[[xf_numfmt[idx]]]
            if (!is.null(code)) codes <- c(codes, code)
          }
        }
        checked <- checked + 1
      }
      if (checked >= n_check + 1) break  # +1 to skip a possible header row
    }
    codes
  }, error = function(e) character(0))
  result
}

#' Robustly read a plate-reader Time column as elapsed seconds
#'
#' Tries, in order: (1) inspect the workbook's actual XML formatting for
#' the time column -- if it's a genuine Excel time/duration format, the raw
#' value openxlsx returns is a fraction of a day and gets multiplied by
#' 86400; if it's a plain number format, the raw value is assumed to
#' already be seconds. (2) If the workbook can't be introspected this way
#' (unusual file, parse failure), falls back to a magnitude heuristic:
#' raw values with a median step under 10 are almost certainly
#' fraction-of-day (a 30-day experiment is still < 30, while a real
#' seconds-denominated sampling interval is essentially never under 10
#' seconds for OD growth curves).
#'
#' @param time_raw Numeric vector as returned by `openxlsx::read.xlsx()`
#'   for the Time column.
#' @param file_path,sheet,col_letter Passed to [susi_excel_cell_formats()]
#'   for the primary detection method.
#' @return List: `time_sec` (numeric, elapsed seconds) and `method`
#'   (`"xml_time_format"`, `"xml_plain_format"`, `"magnitude_fallback"`).
#' @export
susi_resolve_time_seconds <- function(time_raw, file_path = NULL, sheet = NULL, col_letter = "A") {
  codes <- character(0)
  if (!is.null(file_path) && !is.null(sheet) && requireNamespace("xml2", quietly = TRUE)) {
    codes <- susi_excel_cell_formats(file_path, sheet, col_letter)
  }
  if (length(codes) > 0) {
    is_time <- any(vapply(codes, .is_time_format_code, logical(1)))
    if (is_time) {
      return(list(time_sec = time_raw * 86400, method = "xml_time_format"))
    }
    return(list(time_sec = time_raw, method = "xml_plain_format"))
  }
  # Fallback: magnitude heuristic
  step <- suppressWarnings(stats::median(diff(sort(unique(time_raw))), na.rm = TRUE))
  if (!is.na(step) && step > 0 && step < 10) {
    return(list(time_sec = time_raw * 86400, method = "magnitude_fallback"))
  }
  list(time_sec = time_raw, method = "magnitude_fallback")
}
