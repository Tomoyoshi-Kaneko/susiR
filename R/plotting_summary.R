
# ============================================================================
# Additional views, added after user feedback that a single "everything at
# a glance" overlay plot, and bar/point charts of the indices themselves
# (not just a results table), were part of the original workflow.
#
# Neither of these needed anything condition-label-specific: the x-axis is
# just "whatever labels this dataset's conditions are" (categorical), with
# a numeric-MOI-aware upgrade (log-scaled x-axis, viridis gradient colour)
# used only when susi_check_moi_applicable() says the labels support it. So
# this works the same way whether conditions are an MOI series or, say,
# several phage names at a fixed MOI.
# ============================================================================

#' Combined overlay: every condition's mean curve on a single plot
#'
#' A single-panel "everything at a glance" view -- the control curve plus
#' every condition's treated curve, colour-coded, each shown as a mean +/-
#' SD band **across biological replicates** (not across all wells pooled,
#' which would mix well-to-well and day-to-day variability together). If a
#' condition has only one biological replicate, its band collapses to the
#' line (SD = 0 -- nothing to show). If the condition labels parse as a
#' numeric MOI series ([susi_check_moi_applicable()]), colour follows a
#' continuous MOI gradient (and the legend is ordered accordingly);
#' otherwise a plain categorical colour scale is used. No t0/ti/tc markers
#' are drawn here (with many curves at once they would be unreadable) --
#' use [plot_condition()] or [plot_all_conditions()] for that level of
#' per-condition detail.
#'
#' @inheritParams run_susi
#' @param show_error_band If `TRUE` (default), shades +/- 1 SD across
#'   biological replicates around each mean curve.
#' @param condition_colors Optional named character vector mapping condition
#'   label -> colour (a common name like `"red"` or `"steelblue"` works
#'   directly; a hex code like `"#1b9e77"` works too, for anything more
#'   specific), for the categorical (non-MOI) case. Unmatched
#'   conditions fall back to `ggplot2`'s default palette. Ignored when the
#'   conditions form an MOI series (a continuous viridis gradient is used
#'   there by design).
#' @return A `ggplot` object.
#' @export
plot_combined_curves <- function(file_path,
                                  od_sheet = "R", map_sheet = "name",
                                  well_col = "num", condition_col = "sample", bio_rep_col = "bio_rep",
                                  host_col = NULL, host = NULL,
                                  tech_reps_per_bio_rep = 12,
                                  control_label = NULL, conditions = NULL, exclude = NULL, condition_order = NULL,
                                  time_limit_hours = NULL, show_error_band = TRUE, condition_colors = NULL) {
  data <- read_plate_data(file_path, od_sheet, map_sheet, well_col, condition_col, bio_rep_col, tech_reps_per_bio_rep)
  data <- susi_filter_by_host(data, host_col, host, well_col)
  cond_info <- susi_resolve_conditions(data$od_long$condition, control_label, conditions, exclude, condition_order)
  control <- cond_info$control
  conds <- cond_info$conditions

  od_long <- data$od_long
  if (!is.null(time_limit_hours)) od_long <- od_long[od_long$time_h <= time_limit_hours, ]

  ## Mean +/- SD *across biological replicates* at each time point (i.e.
  ## first collapse technical replicates within each bio_rep, then treat
  ## each bio_rep's curve as one observation for the error band).
  per_condition_stats <- function(cond) {
    sub <- od_long[od_long$condition == cond, ]
    per_rep <- stats::aggregate(od ~ time_h + bio_rep, data = sub, FUN = mean, na.rm = TRUE)
    agg <- stats::aggregate(od ~ time_h, data = per_rep,
                             FUN = function(x) c(mean = mean(x), sd = if (length(x) > 1) stats::sd(x) else 0))
    out <- data.frame(time_h = agg$time_h, mean_od = agg$od[, "mean"], sd_od = agg$od[, "sd"])
    out$sd_od[is.na(out$sd_od)] <- 0
    cbind(out, condition = cond)
  }

  df <- do.call(rbind, lapply(c(control, conds), per_condition_stats))
  df$is_control <- df$condition == control
  df$condition <- factor(df$condition, levels = c(control, conds))

  moi_check <- susi_check_moi_applicable(conds)
  p <- ggplot2::ggplot()
  trt <- df[!df$is_control, ]

  if (moi_check$applicable) {
    trt$moi <- moi_check$moi_numeric[as.character(trt$condition)]
    if (show_error_band) {
      p <- p + ggplot2::geom_ribbon(
        data = trt, ggplot2::aes(.data$time_h, ymin = .data$mean_od - .data$sd_od, ymax = .data$mean_od + .data$sd_od,
                                  fill = log10(.data$moi), group = .data$condition), alpha = 0.25, color = NA) +
        ggplot2::scale_fill_viridis_c(name = "log10(MOI)", option = "plasma", guide = "none")
    }
    p <- p +
      ggplot2::geom_line(data = trt, ggplot2::aes(.data$time_h, .data$mean_od, color = log10(.data$moi), group = .data$condition), linewidth = 0.8) +
      ggplot2::scale_color_viridis_c(name = "log10(MOI)", option = "plasma")
  } else {
    if (show_error_band) {
      p <- p + ggplot2::geom_ribbon(
        data = trt, ggplot2::aes(.data$time_h, ymin = .data$mean_od - .data$sd_od, ymax = .data$mean_od + .data$sd_od,
                                  fill = .data$condition, group = .data$condition), alpha = 0.25, color = NA) +
        ggplot2::guides(fill = "none")
      if (!is.null(condition_colors)) p <- p + ggplot2::scale_fill_manual(values = condition_colors)
    }
    p <- p +
      ggplot2::geom_line(data = trt, ggplot2::aes(.data$time_h, .data$mean_od, color = .data$condition), linewidth = 0.8) +
      ggplot2::labs(color = "Condition")
    if (!is.null(condition_colors)) p <- p + ggplot2::scale_color_manual(values = condition_colors)
  }

  ctrl <- df[df$is_control, ]
  if (show_error_band) {
    p <- p + ggplot2::geom_ribbon(data = ctrl, ggplot2::aes(.data$time_h, ymin = .data$mean_od - .data$sd_od, ymax = .data$mean_od + .data$sd_od),
                                   fill = "black", alpha = 0.15)
  }
  p +
    ggplot2::geom_line(data = ctrl, ggplot2::aes(.data$time_h, .data$mean_od),
                        color = "black", linewidth = 1.1, linetype = "solid") +
    ggplot2::annotate("text", x = max(df$time_h), y = ctrl$mean_od[which.max(ctrl$time_h)],
                       label = paste0("  ", control), hjust = 0, size = 3, fontface = "bold") +
    ggplot2::labs(title = "All conditions vs. control (mean \u00b1 SD across biological replicates)",
                  x = "Time (h)", y = "OD") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 40, 5.5, 5.5))
}


#' Colour palette for biological replicates
#'
#' First three replicates get red/green/blue (matching the color coding
#' used in the companion manuscript's figures); further replicates extend
#' with additional distinguishable hues.
#'
#' @param bio_reps Character or numeric vector of biological-replicate IDs.
#' @return Named character vector of hex colours, one per unique ID.
#' @keywords internal
susi_bio_rep_palette <- function(bio_reps) {
  ids <- sort(unique(as.integer(bio_reps)))
  base <- c("#e74c3c", "#27ae60", "#2980b9", "#f39c12", "#8e44ad", "#16a085")
  cols <- if (length(ids) <= length(base)) base[seq_along(ids)] else
    grDevices::colorRampPalette(base)(length(ids))
  stats::setNames(cols, as.character(ids))
}

#' Well-, replicate- and condition-level view of a metric across conditions
#'
#' A single view combining all three levels of the analysis, matching the
#' figure style used in the companion manuscript: every individual well's
#' value as a small circle, each biological replicate's mean as a larger
#' triangle, and the overall condition mean as a black diamond with an
#' error bar (+/- SD across biological replicates). Wells and replicate
#' means share the same red/green/blue-style colour coding by replicate, so
#' within-replicate (well-to-well) and between-replicate (day-to-day)
#' variability are both visible at once, not just collapsed into a single
#' error bar.
#'
#' @inheritParams run_susi
#' @param metrics Which metric(s) to plot. A single metric returns one
#'   `ggplot`; more than one returns a `patchwork` grid (default: all four).
#' @return A `ggplot` object (one metric) or `patchwork` object (several).
#' @export
plot_metric_superplot <- function(file_path, metrics = c("SusI", "VI_local", "SupI", "ti_tc"),
                                   od_sheet = "R", map_sheet = "name",
                                   well_col = "num", condition_col = "sample", bio_rep_col = "bio_rep",
                                   host_col = NULL, host = NULL,
                                   tech_reps_per_bio_rep = 12,
                                   control_label = NULL, conditions = NULL, exclude = NULL, condition_order = NULL,
                                   time_limit_hours = NULL, supi_window_hours = NULL,
                                   params = susi_default_params()) {
  common <- list(od_sheet = od_sheet, map_sheet = map_sheet, well_col = well_col, condition_col = condition_col,
                  bio_rep_col = bio_rep_col, host_col = host_col, tech_reps_per_bio_rep = tech_reps_per_bio_rep,
                  control_label = control_label, conditions = conditions, exclude = exclude, condition_order = condition_order,
                  time_limit_hours = time_limit_hours, supi_window_hours = supi_window_hours, params = params,
                  verbose = FALSE)

  wells_res <- do.call(run_susi, c(list(file_path = file_path, calculation_method = "individual_wells"), common))
  reps_res  <- do.call(run_susi, c(list(file_path = file_path, calculation_method = "biological_replicates"), common))

  if (!is.null(host_col)) {
    if (is.null(host)) stop("`host_col` was given but `host` was not. Available hosts: ",
                             paste(names(wells_res$conditions), collapse = ", "), call. = FALSE)
    level_order <- wells_res$conditions[[host]]
    wells_summary <- wells_res$summary[wells_res$summary$host == host, ]
    reps_summary  <- reps_res$summary[reps_res$summary$host == host, ]
  } else {
    level_order <- wells_res$conditions
    wells_summary <- wells_res$summary
    reps_summary <- reps_res$summary
  }

  wells_df <- wells_summary
  wells_df$condition <- factor(wells_df$condition, levels = level_order)
  wells_df$bio_rep <- factor(wells_df$bio_rep)

  reps_df <- reps_summary[!is.na(reps_summary$bio_rep), ]
  reps_df$condition <- factor(reps_df$condition, levels = level_order)
  reps_df$bio_rep <- factor(reps_df$bio_rep)

  cond_df <- reps_summary[is.na(reps_summary$bio_rep), ]
  cond_df$condition <- factor(cond_df$condition, levels = level_order)

  pal <- susi_bio_rep_palette(levels(wells_df$bio_rep))
  dodge_w <- 0.6
  jd <- ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_w, seed = 1)
  dg <- ggplot2::position_dodge(width = dodge_w)

  one_panel <- function(metric) {
    se_col <- paste0(metric, "_se")
    cond_df$se <- if (se_col %in% names(cond_df)) cond_df[[se_col]] else NA_real_

    ggplot2::ggplot() +
      ggplot2::geom_point(data = wells_df, ggplot2::aes(x = .data$condition, y = .data[[metric]], color = .data$bio_rep),
                           position = jd, shape = 16, size = 1.7, alpha = 0.55) +
      ggplot2::geom_point(data = reps_df, ggplot2::aes(x = .data$condition, y = .data[[metric]], color = .data$bio_rep),
                           position = dg, shape = 17, size = 3.2, stroke = 0.6) +
      ggplot2::geom_errorbar(data = cond_df, ggplot2::aes(x = .data$condition, ymin = .data[[metric]] - .data$se, ymax = .data[[metric]] + .data$se),
                              width = 0.18, linewidth = 0.7, color = "black") +
      ggplot2::geom_point(data = cond_df, ggplot2::aes(x = .data$condition, y = .data[[metric]]),
                           shape = 23, size = 3.2, fill = "black", color = "black") +
      ggplot2::scale_color_manual(values = pal, name = "Bio rep") +
      ggplot2::labs(title = metric, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                      plot.title = ggplot2::element_text(face = "bold"))
  }

  panels <- lapply(metrics, one_panel)
  if (length(panels) == 1) return(panels[[1]] + ggplot2::labs(caption = "Circles = individual wells; triangles = biological-replicate means; black diamond \u00b1 error bar = condition mean \u00b1 SE"))
  patchwork::wrap_plots(panels, ncol = min(2, length(panels)), guides = "collect") +
    patchwork::plot_annotation(caption = "Circles = individual wells; triangles = biological-replicate means; black diamond \u00b1 error bar = condition mean \u00b1 SE")
}



#' Bar/point charts of SusI, VI, SupI and ti/tc across conditions
#'
#' Turns the numbers in `run_susi()$summary` into one small chart per
#' metric (2x2 grid by default). Works the same way regardless of what the
#' condition labels mean: if [susi_check_moi_applicable()] recognizes them
#' as an MOI dilution series, points are plotted against log10(MOI) with a
#' connecting line (dose-response style); otherwise each condition is
#' simply a labelled bar/point in the order [run_susi()] analyzed them.
#' When `calculation_method = "biological_replicates"` was used, error bars
#' show +/- 1 SE (standard error) across biological replicates, matching
#' the companion manuscript's figures. For a view that also shows the
#' underlying wells and per-replicate means, see [plot_metric_superplot()].
#'
#' @param result The list returned by [run_susi()].
#' @param metrics Which metrics to plot. Default all four.
#' @param condition_colors Optional named character vector mapping condition
#'   label -> colour (any value `ggplot2` understands, e.g. `"#c0392b"` or
#'   `"steelblue"`), for the categorical (non-MOI) bar-chart case. Unnamed
#'   / unmatched conditions fall back to the default colour. Ignored when
#'   the conditions form an MOI dose-response series (a single line/point
#'   colour is used there, matching the plot's continuous-fit style).
#' @return A `ggplot` object (single metric) or `patchwork` object (multiple).
#' @export
plot_metrics_summary <- function(result, metrics = c("SusI", "VI_local", "SupI", "ti_tc"), condition_colors = NULL) {
  summ <- result$summary
  cond_level <- if ("well" %in% names(summ)) {
    dplyr::summarise(dplyr::group_by(summ, .data$condition),
                      dplyr::across(dplyr::all_of(metrics), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
                      dplyr::across(dplyr::all_of(metrics), ~ {
                        x <- .x[!is.na(.x)]
                        if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else 0
                      }, .names = "{.col}_se"),
                      .groups = "drop")
  } else if ("bio_rep" %in% names(summ) && any(!is.na(summ$bio_rep))) {
    summ[is.na(summ$bio_rep), ]
  } else {
    summ
  }
  cond_level$condition <- factor(cond_level$condition, levels = result$conditions)

  moi_check <- susi_check_moi_applicable(result$conditions)
  dose_response <- moi_check$applicable

  one_panel <- function(metric) {
    se_col <- paste0(metric, "_se")
    has_se <- se_col %in% names(cond_level) && any(!is.na(cond_level[[se_col]]))
    d <- cond_level
    d$value <- d[[metric]]
    d$se <- if (has_se) d[[se_col]] else NA_real_

    if (dose_response) {
      d$moi <- moi_check$moi_numeric[as.character(d$condition)]
      p <- ggplot2::ggplot(d, ggplot2::aes(x = log10(.data$moi), y = .data$value)) +
        ggplot2::geom_line(color = "#2c3e50") +
        ggplot2::geom_point(size = 2.2, color = "#c0392b") +
        ggplot2::labs(x = "log10(MOI)")
    } else {
      p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$condition, y = .data$value))
      if (is.null(condition_colors)) {
        p <- p + ggplot2::geom_col(fill = "#c0392b", width = 0.6, alpha = 0.85)
      } else {
        p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data$condition), width = 0.6, alpha = 0.9) +
          ggplot2::scale_fill_manual(values = condition_colors, na.value = "#c0392b") +
          ggplot2::guides(fill = "none")
      }
      p <- p + ggplot2::labs(x = NULL) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    }
    if (has_se) {
      p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$value - .data$se, ymax = .data$value + .data$se), width = 0.15)
    }
    p + ggplot2::labs(title = metric, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }

  panels <- lapply(metrics, one_panel)
  if (length(panels) == 1) return(panels[[1]])
  patchwork::wrap_plots(panels, ncol = min(2, length(panels)))
}
