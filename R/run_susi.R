# ============================================================================
# run_susi.R -- top-level entry point
# ============================================================================

#' Compute SusI, VI, SupI and ti/tc for every condition in a plate dataset
#'
#' This is the main entry point. It reads the workbook, automatically
#' resolves which wells belong to which condition and biological replicate
#' (see [read_plate_data()] / [susi_resolve_bio_rep()]), detects t0/ti/tc,
#' and computes all four indices -- with **no assumption about the number
#' of conditions, their labels, or the number of biological replicates**.
#'
#' @inheritParams read_plate_data
#' @param control_label,conditions,exclude See [susi_resolve_conditions()].
#' @param calculation_method One of `"biological_replicates"` (default;
#'   mean curve per biological replicate, then averaged -- matches the
#'   validated pipeline described in the companion manuscript),
#'   `"individual_wells"` (every well analyzed independently against the
#'   pooled control mean; no averaging), or `"overall_mean"` (a single mean
#'   curve per condition).
#' @param supi_window_hours Upper bound (hours) of the SupI integration
#'   window. Defaults to `time_limit_hours` if set (matching how the
#'   companion manuscript's actual pipeline was run -- `calculate_suppression_index_v20()`'s
#'   own hard-coded default of 30 h was never used in practice, since
#'   `process_phage_data_v20()` always passed its own `time_limit_hours`
#'   through instead), else `30`.
#' @param time_limit_hours Optional overall upper time limit (hours) applied
#'   throughout (detection and integration). `NULL` = use the full time course.
#' @param params Detection parameters, see [susi_default_params()].
#' @return A list:
#'   \item{summary}{One row per condition (plus one row per condition x
#'     bio_rep when `calculation_method = "biological_replicates"`) with
#'     `SusI`, `VI_local`, `SupI`, `ti_tc`, `t0`, `ti`, `tc`.}
#'   \item{global_vi}{`VI` and `MV50` across the dilution series, or `NA`
#'     with an explanatory `note` if the conditions do not form an MOI
#'     series (see [susi_check_moi_applicable()]).}
#'   \item{diagnostics}{Per-well/per-curve detection outcome table.}
#'   \item{control}{The control condition label used.}
#'   \item{params}{The detection parameters used.}
#' @export
run_susi <- function(file_path,
                      od_sheet = "R", map_sheet = "name",
                      well_col = "num", condition_col = "sample", bio_rep_col = "bio_rep",
                      tech_reps_per_bio_rep = 12,
                      control_label = NULL, conditions = NULL, exclude = NULL,
                      calculation_method = c("biological_replicates", "individual_wells", "overall_mean"),
                      time_limit_hours = NULL,
                      supi_window_hours = NULL,
                      params = susi_default_params(),
                      verbose = TRUE) {
  calculation_method <- match.arg(calculation_method)
  if (!is.null(time_limit_hours)) params$time_limit_hours <- time_limit_hours
  if (is.null(supi_window_hours)) supi_window_hours <- if (!is.null(params$time_limit_hours)) params$time_limit_hours else 30

  data <- read_plate_data(file_path, od_sheet, map_sheet, well_col, condition_col,
                           bio_rep_col, tech_reps_per_bio_rep)
  cond_info <- susi_resolve_conditions(data$mapping[[condition_col]], control_label, conditions, exclude)
  control <- cond_info$control
  conds <- cond_info$conditions
  if (verbose) message("Conditions to analyze (", length(conds), "): ", paste(conds, collapse = ", "))

  od_long <- data$od_long
  time_h <- data$time_h

  wide_curve <- function(df, value_col = "od") {
    stats::aggregate(as.formula(paste(value_col, "~ time_h")), data = df, FUN = mean, na.rm = TRUE)
  }

  control_all <- od_long[od_long$condition == control, ]
  control_mean_overall <- wide_curve(control_all)
  tc_overall <- susi_detect_tc(control_mean_overall$time_h, control_mean_overall$od, params)
  ## Effective upper integration bound used throughout (SusI denominator, ti
  ## fallback, VI window): min(tc, time_limit_hours), matching the legacy
  ## scripts' `effective_end_time`. tc itself is always a single,
  ## dataset-wide value detected once from the pooled control curve --
  ## never per-condition or per-biological-replicate.
  effective_end_time <- if (!is.null(params$time_limit_hours)) min(tc_overall$tc, params$time_limit_hours) else tc_overall$tc

  summary_rows <- list()
  diag_rows <- list()

  add_diag <- function(condition, bio_rep, well, t0_res, ti_res, susi_val) {
    diag_rows[[length(diag_rows) + 1]] <<- data.frame(
      condition = condition, bio_rep = bio_rep, well = if (is.null(well)) NA_character_ else well,
      t0 = t0_res$t0, t0_method = t0_res$method,
      ti = ti_res$ti, ti_detected = ti_res$detected,
      status = dplyr::case_when(
        is.na(t0_res$t0) ~ "t0_detection_failed",
        is.na(susi_val) ~ "susi_calculation_failed",
        TRUE ~ "success"
      ),
      stringsAsFactors = FALSE
    )
  }

  if (calculation_method == "overall_mean") {
    for (cond in conds) {
      treated_mean <- wide_curve(od_long[od_long$condition == cond, ])
      t0_res <- susi_detect_t0(treated_mean$time_h, treated_mean$od, params)
      ti_res <- if (is.na(t0_res$t0)) list(ti = NA_real_, detected = NA) else
        susi_detect_ti(treated_mean$time_h, treated_mean$od, t0_res$t0, params, fallback_limit = effective_end_time)

      c_al <- control_mean_overall$od[match(time_h, control_mean_overall$time_h)]
      p_al <- treated_mean$od[match(time_h, treated_mean$time_h)]

      susi_val <- susi_calc_susi(time_h, c_al, p_al, t0_res$t0, ti_res$ti, effective_end_time)
      vi_val <- susi_calc_vi_local(time_h, c_al, p_al, effective_end_time)
      supi_val <- susi_calc_supi(time_h, c_al, p_al, supi_window_hours)
      tt_val <- susi_calc_time_ratio(ti_res$ti, tc_overall$tc)

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        condition = cond, bio_rep = NA_integer_,
        SusI = susi_val, VI_local = vi_val, SupI = supi_val, ti_tc = tt_val,
        t0 = t0_res$t0, ti = ti_res$ti, tc = tc_overall$tc,
        n_wells = length(unique(od_long$well[od_long$condition == cond])),
        stringsAsFactors = FALSE
      )
      add_diag(cond, NA_integer_, NULL, t0_res, ti_res, susi_val)
    }

  } else if (calculation_method == "biological_replicates") {
    for (cond in conds) {
      bio_reps <- sort(unique(od_long$bio_rep[od_long$condition == cond]))
      per_rep <- list()
      for (br in bio_reps) {
        treated_mean <- wide_curve(od_long[od_long$condition == cond & od_long$bio_rep == br, ])
        control_mean_br <- wide_curve(control_all[control_all$bio_rep == br, ])
        if (nrow(control_mean_br) == 0) {
          if (verbose) message("  [", cond, " / bio_rep ", br, "] no matching control replicate -- skipped.")
          next
        }
        t0_res <- susi_detect_t0(treated_mean$time_h, treated_mean$od, params)
        ti_res <- if (is.na(t0_res$t0)) list(ti = NA_real_, detected = NA) else
          susi_detect_ti(treated_mean$time_h, treated_mean$od, t0_res$t0, params, fallback_limit = effective_end_time)

        # SusI numerator (t0->ti area) uses this bio-rep's own control curve;
        # the denominator (0->effective_end_time) uses the POOLED control
        # curve across all replicates, matching the legacy scripts exactly
        # (a single, shared normalization constant, not one per replicate).
        c_al_bio <- control_mean_br$od[match(time_h, control_mean_br$time_h)]
        c_al_pooled <- control_mean_overall$od[match(time_h, control_mean_overall$time_h)]
        p_al <- treated_mean$od[match(time_h, treated_mean$time_h)]

        susi_val <- susi_calc_susi(time_h, c_al_bio, p_al, t0_res$t0, ti_res$ti, effective_end_time,
                                    C_denom = c_al_pooled)
        vi_val   <- susi_calc_vi_local(time_h, c_al_bio, p_al, effective_end_time)
        supi_val <- susi_calc_supi(time_h, c_al_bio, p_al, supi_window_hours)
        tt_val   <- susi_calc_time_ratio(ti_res$ti, tc_overall$tc)

        n_w <- length(unique(od_long$well[od_long$condition == cond & od_long$bio_rep == br]))
        per_rep[[length(per_rep) + 1]] <- data.frame(
          condition = cond, bio_rep = br, SusI = susi_val, VI_local = vi_val, SupI = supi_val, ti_tc = tt_val,
          t0 = t0_res$t0, ti = ti_res$ti, tc = tc_overall$tc, n_wells = n_w, stringsAsFactors = FALSE
        )
        add_diag(cond, br, NULL, t0_res, ti_res, susi_val)
      }
      if (length(per_rep) == 0) next
      rep_df <- dplyr::bind_rows(per_rep)
      summary_rows[[length(summary_rows) + 1]] <- rep_df
      safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      safe_se <- function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else 0
      }
      agg <- data.frame(
        condition = cond, bio_rep = NA_integer_,
        SusI = safe_mean(rep_df$SusI), VI_local = safe_mean(rep_df$VI_local),
        SupI = safe_mean(rep_df$SupI), ti_tc = safe_mean(rep_df$ti_tc),
        t0 = safe_mean(rep_df$t0), ti = safe_mean(rep_df$ti), tc = safe_mean(rep_df$tc),
        n_wells = sum(rep_df$n_wells),
        # Standard ERROR across biological replicates (sd / sqrt(n)) -- this
        # is what the companion manuscript's figures plot as error bars, not
        # raw SD (which is ~1.7x larger at n=3 and was the main reason this
        # package's error bars looked inflated versus the paper's figures).
        SusI_se = safe_se(rep_df$SusI),
        VI_local_se = safe_se(rep_df$VI_local),
        SupI_se = safe_se(rep_df$SupI),
        ti_tc_se = safe_se(rep_df$ti_tc),
        n_bio_reps = nrow(rep_df),
        stringsAsFactors = FALSE
      )
      summary_rows[[length(summary_rows) + 1]] <- agg
    }

  } else { # individual_wells
    for (cond in conds) {
      wells <- unique(od_long$well[od_long$condition == cond])
      for (w in wells) {
        well_curve <- od_long[od_long$well == w, c("time_h", "od")]
        well_curve <- well_curve[match(time_h, well_curve$time_h), ]
        t0_res <- susi_detect_t0(well_curve$time_h, well_curve$od, params)
        ti_res <- if (is.na(t0_res$t0)) list(ti = NA_real_, detected = NA) else
          susi_detect_ti(well_curve$time_h, well_curve$od, t0_res$t0, params, fallback_limit = effective_end_time)

        c_al <- control_mean_overall$od[match(time_h, control_mean_overall$time_h)]
        p_al <- well_curve$od

        susi_val <- susi_calc_susi(time_h, c_al, p_al, t0_res$t0, ti_res$ti, effective_end_time)
        vi_val   <- susi_calc_vi_local(time_h, c_al, p_al, effective_end_time)
        supi_val <- susi_calc_supi(time_h, c_al, p_al, supi_window_hours)
        tt_val   <- susi_calc_time_ratio(ti_res$ti, tc_overall$tc)

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          condition = cond, bio_rep = od_long$bio_rep[od_long$well == w][1], well = w,
          SusI = susi_val, VI_local = vi_val, SupI = supi_val, ti_tc = tt_val,
          t0 = t0_res$t0, ti = ti_res$ti, tc = tc_overall$tc, n_wells = 1, stringsAsFactors = FALSE
        )
        add_diag(cond, od_long$bio_rep[od_long$well == w][1], w, t0_res, ti_res, susi_val)
      }
    }
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  diag_df <- if (length(diag_rows) > 0) dplyr::bind_rows(diag_rows) else data.frame()

  # One row of local VI per *condition* (averaging over wells for
  # individual_wells; the aggregate row already exists for the other two
  # methods) is what the global MOI dose-response fit needs.
  cond_level_vi <- switch(calculation_method,
    overall_mean          = stats::setNames(summary_df$VI_local, summary_df$condition),
    biological_replicates = { agg <- summary_df[is.na(summary_df$bio_rep), ]; stats::setNames(agg$VI_local, agg$condition) },
    individual_wells       = { m <- stats::aggregate(VI_local ~ condition, data = summary_df, FUN = mean, na.rm = TRUE); stats::setNames(m$VI_local, m$condition) }
  )

  moi_check <- susi_check_moi_applicable(conds)
  global_vi <- list(VI = NA_real_, MV50 = NA_real_, applicable = FALSE,
                     note = "Conditions do not form a numeric MOI dilution series (need >=3 distinct numeric MOI labels); global VI/MV50 skipped. Local, per-condition VI is still reported in `summary`.")
  if (moi_check$applicable) {
    moi_vec <- moi_check$moi_numeric[names(cond_level_vi)]
    fit <- susi_calc_global_vi(moi_vec, cond_level_vi)
    global_vi <- list(VI = fit$VI, MV50 = fit$MV50, applicable = TRUE, note = NULL)
  }

  list(summary = summary_df, global_vi = global_vi, diagnostics = diag_df,
       control = control, conditions = conds, calculation_method = calculation_method, params = params)
}
