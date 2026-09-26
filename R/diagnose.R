# ============================================================================
# diagnose.R -- rapid, metric-free detection diagnostics for parameter
#               tuning (per-well t0/ti detection outcome), independent of
#               how many conditions/replicates are present.
# ============================================================================

#' Diagnose t0/ti detection success across every well, without computing indices
#'
#' Much cheaper than a full [run_susi()] call: runs only the detection step
#' (§2.2 of the Application Note) on every individual well and classifies
#' the outcome, so detection parameters can be tuned quickly before running
#' the full analysis. Useful the same way regardless of how many conditions
#' or replicates the dataset has.
#'
#' @inheritParams run_susi
#' @param output_excel Optional path; if given, writes an annotated Excel
#'   workbook (one row per well, one sheet per condition) alongside a
#'   Summary sheet.
#' @return A data frame, one row per well: `condition`, `bio_rep`, `well`,
#'   `t0`, `t0_method`, `ti`, `ti_detected`, `status` (plus `host` when
#'   `host_col` is used).
#' @export
diagnose_conditions <- function(file_path,
                                 od_sheet = "R", map_sheet = "name",
                                 well_col = "well", condition_col = "sample", bio_rep_col = "bio_rep",
                                 host_col = NULL,
                                 tech_reps_per_bio_rep = 12,
                                 control_label = NULL, conditions = NULL, exclude = NULL, condition_order = NULL,
                                 params = susi_default_params(),
                                 output_excel = NULL,
                                 verbose = TRUE) {
  data <- read_plate_data(file_path, od_sheet, map_sheet, well_col, condition_col,
                           bio_rep_col, tech_reps_per_bio_rep)
  od_long <- data$od_long
  time_h <- data$time_h

  use_host <- !is.null(host_col) && host_col %in% names(data$mapping)
  if (use_host) {
    host_map <- stats::setNames(data$mapping[[host_col]], data$mapping[[well_col]])
    od_long$host <- host_map[od_long$well]
    od_long <- od_long[!is.na(od_long$host) & nzchar(trimws(od_long$host)), ]
    all_conds <- unique(od_long$condition)
    if (!is.null(exclude)) all_conds <- setdiff(all_conds, exclude)
  } else {
    cond_info <- susi_resolve_conditions(od_long$condition, control_label, conditions, exclude, condition_order)
    all_conds <- c(cond_info$control, cond_info$conditions)
  }

  rows <- list()

  for (cond in all_conds) {
    wells <- unique(od_long$well[od_long$condition == cond])
    for (w in wells) {
      wc <- od_long[od_long$well == w, c("time_h", "od")]
      wc <- wc[match(time_h, wc$time_h), ]
      t0_res <- susi_detect_t0(wc$time_h, wc$od, params)
      ti_res <- if (is.na(t0_res$t0)) list(ti = NA_real_, detected = NA) else
        susi_detect_ti(wc$time_h, wc$od, t0_res$t0, params)
      status <- if (is.na(t0_res$t0)) "t0_detection_failed" else if (isFALSE(ti_res$detected)) "ti_not_detected" else "success"
      row <- data.frame(
        condition = cond, bio_rep = od_long$bio_rep[od_long$well == w][1], well = w,
        t0 = t0_res$t0, t0_method = t0_res$method,
        ti = ti_res$ti, ti_detected = ti_res$detected,
        status = status, stringsAsFactors = FALSE
      )
      if (use_host) row$host <- od_long$host[od_long$well == w][1]
      rows[[length(rows) + 1]] <- row
    }
  }
  out <- dplyr::bind_rows(rows)
  if (use_host && nrow(out) > 0) out <- out[, c("host", setdiff(names(out), "host")), drop = FALSE]

  if (verbose) {
    tab <- table(out$condition, out$status)
    message("Detection summary (rows = condition, cols = status):")
    print(tab)
  }

  if (!is.null(output_excel)) {
    wb <- openxlsx::createWorkbook()
    summary_tab <- stats::aggregate(well ~ condition + status, data = out, FUN = length)
    names(summary_tab)[3] <- "n_wells"
    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_tab)
    for (cond in all_conds) {
      sheet <- substr(gsub("[^A-Za-z0-9_]", "_", cond), 1, 31)
      openxlsx::addWorksheet(wb, sheet)
      openxlsx::writeData(wb, sheet, out[out$condition == cond, ])
    }
    non_success <- out[out$status != "success", ]
    if (nrow(non_success) > 0) {
      openxlsx::addWorksheet(wb, "Non_Success")
      openxlsx::writeData(wb, "Non_Success", non_success)
    }
    openxlsx::saveWorkbook(wb, output_excel, overwrite = TRUE)
    if (verbose) message("Diagnostic workbook written to ", output_excel)
  }

  out
}
