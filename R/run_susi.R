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
#' If a dataset contains **more than one host** (e.g. several bacterial
#' strains sharing a plate, each with its own phage-free control), pass
#' `host_col` naming a column in the mapping sheet that labels every well
#' (control and treated alike) with which host it belongs to. Each host is
#' then analyzed independently -- its own control resolved, its own
#' conditions computed, its own global MOI fit attempted -- and the results
#' are combined with an added `host` column. Without `host_col` (the
#' default), behavior is exactly as before: one shared control for the
#' whole dataset.
#'
#' @inheritParams read_plate_data
#' @param control_label,conditions,exclude See [susi_resolve_conditions()].
#'   When `host_col` is used, `control_label` may also be a named
#'   character vector (one entry per host, e.g.
#'   `c(A5940 = "Host A5940", VISA = "Host VISA")`) if the control isn't
#'   auto-detectable within every host's own subset; a single string is
#'   still fine if every host's control shares a common recognizable label
#'   in [susi_resolve_conditions()]'s auto-detection.
#' @param host_col Optional: name of a column in the mapping sheet that
#'   assigns every well to a host/strain group. When given, susiR analyzes
#'   each host's wells separately (its own control, its own conditions),
#'   rather than assuming the whole sheet shares one control. Wells with a
#'   blank/`NA` host value are dropped (use this for shared blanks like a
#'   media-only control that isn't tied to one host, and exclude them via
#'   `exclude` within each host's own conditions as needed).
#' @param calculation_method One of `"biological_replicates"` (default;
#'   mean curve per biological replicate, then averaged -- matches the
#'   validated pipeline described in the companion manuscript),
#'   `"individual_wells"` (every well analyzed independently against the
#'   pooled control mean; no averaging), or `"overall_mean"` (a single mean
#'   curve per condition).
#' @param supi_window_hours Upper bound (hours) of the SupI integration
#'   window. Defaults to `time_limit_hours` if set, else `30`.
#' @param time_limit_hours Optional overall upper time limit (hours) applied
#'   throughout (detection and integration). `NULL` = use the full time course.
#' @param params Detection parameters, see [susi_default_params()].
#' @return A list:
#'   \item{summary}{One row per condition (plus one row per condition x
#'     bio_rep when `calculation_method = "biological_replicates"`) with
#'     `SusI`, `VI_local`, `SupI`, `ti_tc`, `t0`, `ti`, `tc`, and (when
#'     `host_col` is used) `host`.}
#'   \item{global_vi}{`VI` and `MV50` across the dilution series (a single
#'     list, or -- when `host_col` is used -- a named list of one such list
#'     per host), or `NA` with an explanatory `note` if the conditions do
#'     not form an MOI series (see [susi_check_moi_applicable()]).}
#'   \item{diagnostics}{Per-well/per-curve detection outcome table.}
#'   \item{control}{The control condition label used (or, with `host_col`,
#'     a named vector of one per host).}
#'   \item{params}{The detection parameters used.}
#' @export
run_susi <- function(file_path,
                      od_sheet = "R", map_sheet = "name",
                      well_col = "well", condition_col = "sample", bio_rep_col = "bio_rep",
                      host_col = NULL,
                      tech_reps_per_bio_rep = 12,
                      control_label = NULL, conditions = NULL, exclude = NULL, condition_order = NULL,
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

  use_host <- !is.null(host_col) && host_col %in% names(data$mapping)
  if (!is.null(host_col) && !use_host) {
    stop("`host_col` = '", host_col, "' does not match any column in the mapping sheet (",
         paste(names(data$mapping), collapse = ", "), ").", call. = FALSE)
  }

  if (!use_host) {
    core <- .run_susi_core(data$od_long, data$time_h, condition_col,
                            control_label, conditions, exclude, condition_order,
                            calculation_method, supi_window_hours, params, verbose)
    return(list(summary = core$summary, global_vi = core$global_vi, diagnostics = core$diagnostics,
                control = core$control, conditions = core$conditions,
                calculation_method = calculation_method, params = params))
  }

  ## --- multi-host path -----------------------------------------------------
  host_map <- stats::setNames(data$mapping[[host_col]], data$mapping[[well_col]])
  od_long_all <- data$od_long
  od_long_all$host <- host_map[od_long_all$well]

  hosts <- unique(od_long_all$host)
  hosts <- hosts[!is.na(hosts) & nzchar(trimws(as.character(hosts)))]
  if (length(hosts) == 0) {
    stop("`host_col` = '", host_col, "' exists but every well's value is blank/NA.", call. = FALSE)
  }
  ## Sanity check: a common mistake is pointing host_col at the condition
  ## column itself (or another column that varies with condition), which
  ## makes every "host" a group of one, with no real control -- this fails
  ## deep inside the per-host computation with a confusing, generic error
  ## (e.g. "no rows to aggregate"). Catch it here with an explicit message.
  degenerate <- vapply(hosts, function(h) {
    length(unique(od_long_all$condition[!is.na(od_long_all$host) & od_long_all$host == h])) <= 1
  }, logical(1))
  if (any(degenerate)) {
    stop("`host_col` = '", host_col, "' produces group(s) with only one distinct condition label ",
         "(e.g. '", hosts[degenerate][1], "'), so no control-vs-treatment comparison is possible within it. ",
         "This usually means `host_col` was set to the condition/well column itself rather than a ",
         "separate column that groups several conditions (including their shared control) under one host/strain.",
         call. = FALSE)
  }
  if (verbose) message("Hosts found (", length(hosts), "): ", paste(hosts, collapse = ", "))

  control_per_host <- if (is.null(control_label)) {
    stats::setNames(rep(list(NULL), length(hosts)), hosts)
  } else if (!is.null(names(control_label))) {
    stats::setNames(as.list(control_label)[hosts], hosts)
  } else {
    stats::setNames(rep(list(control_label), length(hosts)), hosts)
  }

  all_summary <- list(); all_diag <- list(); all_global_vi <- list(); all_control <- character(0)

  for (h in hosts) {
    if (verbose) message("--- Host: ", h, " ---")
    od_h <- od_long_all[!is.na(od_long_all$host) & od_long_all$host == h, ]
    core <- tryCatch(
      .run_susi_core(od_h, sort(unique(od_h$time_h)), condition_col,
                      control_per_host[[h]], conditions, exclude, condition_order,
                      calculation_method, supi_window_hours, params, verbose),
      error = function(e) {
        warning("Host '", h, "' failed and was skipped: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(core)) next
    core$summary$host <- h
    core$diagnostics$host <- h
    all_summary[[h]] <- core$summary
    all_diag[[h]] <- core$diagnostics
    all_global_vi[[h]] <- core$global_vi
    all_control[h] <- core$control
  }

  summary_df <- dplyr::bind_rows(all_summary)
  diag_df <- dplyr::bind_rows(all_diag)
  if (nrow(summary_df) > 0) summary_df <- summary_df[, c("host", setdiff(names(summary_df), "host")), drop = FALSE]

  list(summary = summary_df, global_vi = all_global_vi, diagnostics = diag_df,
       control = all_control, conditions = lapply(all_summary, function(x) unique(x$condition)),
       calculation_method = calculation_method, params = params)
}

#' @keywords internal
.run_susi_core <- function(od_long, time_h, condition_col,
                            control_label, conditions, exclude, condition_order,
                            calculation_method, supi_window_hours, params, verbose) {
  cond_info <- susi_resolve_conditions(od_long$condition, control_label, conditions, exclude, condition_order)
  control <- cond_info$control
  conds <- cond_info$conditions
  if (verbose) message("Conditions to analyze (", length(conds), "): ", paste(conds, collapse = ", "))

  wide_curve <- function(df, value_col = "od") {
    stats::aggregate(as.formula(paste(value_col, "~ time_h")), data = df, FUN = mean, na.rm = TRUE)
  }

  control_all <- od_long[od_long$condition == control, ]
  control_mean_overall <- wide_curve(control_all)
  tc_overall <- susi_detect_tc(control_mean_overall$time_h, control_mean_overall$od, params)
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
       control = control, conditions = conds)
}
