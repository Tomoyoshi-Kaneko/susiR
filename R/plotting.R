# ============================================================================
# plotting.R -- visual diagnostics: OD curves with detected t0/ti/tc overlaid,
#               and the SusI numerator region shaded, for a single condition
#               or as a one-page grid across every condition in the dataset.
#               Uses the exact same detection functions as run_susi(), so a
#               plot for a given condition/bio_rep always matches the numbers
#               run_susi() reports for that same condition/bio_rep.
# ============================================================================

#' Build one diagnostic panel: control + treated OD curves with t0/ti/tc marked
#'
#' @param time_h Numeric vector of time points (hours), shared by both curves.
#' @param control_od,treated_od Numeric OD vectors, same length as `time_h`.
#' @param t0,ti,tc Detected time points (hours); any may be `NA`.
#' @param title Panel title (e.g. the condition label).
#' @param subtitle Optional panel subtitle (e.g. "bio_rep 2", or the computed SusI value).
#' @param shade_susi If `TRUE` (default), shades the area between the control
#'   and treated curves over `[t0, ti]` -- the SusI numerator.
#' @return A `ggplot` object.
#' @export
susi_plot_curve <- function(time_h, control_od, treated_od, t0 = NA, ti = NA, tc = NA,
                             title = NULL, subtitle = NULL, shade_susi = TRUE) {
  df <- data.frame(
    time_h = rep(time_h, 2),
    od = c(control_od, treated_od),
    curve = rep(c("Control", "Treated"), each = length(time_h))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time_h, y = .data$od, color = .data$curve))

  if (shade_susi && !is.na(t0) && !is.na(ti)) {
    band <- data.frame(time_h = time_h, control_od = control_od, treated_od = treated_od)
    band <- band[band$time_h >= t0 & band$time_h <= ti, ]
    if (nrow(band) >= 2) {
      p <- p + ggplot2::geom_ribbon(
        data = band,
        ggplot2::aes(x = .data$time_h, ymin = .data$treated_od, ymax = .data$control_od),
        inherit.aes = FALSE, fill = "grey70", alpha = 0.45
      )
    }
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c(Control = "#2c3e50", Treated = "#c0392b"))

  mark <- function(p, x, label, linetype, vjust) {
    if (is.na(x)) return(p)
    p +
      ggplot2::geom_vline(xintercept = x, linetype = linetype, color = "grey40", linewidth = 0.5) +
      ggplot2::annotate("text", x = x, y = Inf, label = label, vjust = vjust, hjust = -0.15,
                         size = 3, color = "grey30")
  }
  p <- mark(p, t0, "t0", "dashed", vjust = 1.3)
  p <- mark(p, ti, "ti", "dashed", vjust = 2.6)
  p <- mark(p, tc, "tc", "dotted", vjust = 3.9)

  p +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Time (h)", y = "OD", color = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
}

#' Diagnostic plot for one condition
#'
#' Re-derives the control/treated curves and re-runs detection exactly as
#' [run_susi()] does, for a single condition (and, for
#' `calculation_method = "biological_replicates"`, a single biological
#' replicate if `bio_rep` is given, otherwise every replicate is drawn as
#' separate panels side by side).
#'
#' @inheritParams run_susi
#' @param condition Condition label to plot (must be one of the conditions
#'   resolved by [susi_resolve_conditions()], i.e. not the control).
#' @param bio_rep Optional single biological-replicate ID to restrict to
#'   (only used when `calculation_method = "biological_replicates"`).
#' @return A `ggplot` object (single condition/replicate) or, if multiple
#'   replicates are plotted at once, a combined `patchwork` object.
#' @export
plot_condition <- function(file_path, condition,
                            od_sheet = "R", map_sheet = "name",
                            well_col = "num", condition_col = "sample", bio_rep_col = "bio_rep",
                            tech_reps_per_bio_rep = 12,
                            control_label = NULL, bio_rep = NULL,
                            calculation_method = c("biological_replicates", "overall_mean"),
                            time_limit_hours = NULL,
                            params = susi_default_params()) {
  calculation_method <- match.arg(calculation_method)
  if (!is.null(time_limit_hours)) params$time_limit_hours <- time_limit_hours

  data <- read_plate_data(file_path, od_sheet, map_sheet, well_col, condition_col, bio_rep_col, tech_reps_per_bio_rep)
  cond_info <- susi_resolve_conditions(data$mapping[[condition_col]], control_label)
  control <- cond_info$control
  if (!condition %in% cond_info$conditions) {
    stop("'", condition, "' is not among the resolved conditions: ", paste(cond_info$conditions, collapse = ", "), call. = FALSE)
  }

  od_long <- data$od_long
  time_h <- data$time_h
  wide_curve <- function(df) stats::aggregate(od ~ time_h, data = df, FUN = mean, na.rm = TRUE)
  control_all <- od_long[od_long$condition == control, ]

  if (calculation_method == "overall_mean") {
    ctrl <- wide_curve(control_all)
    trt  <- wide_curve(od_long[od_long$condition == condition, ])
    tc_res <- susi_detect_tc(ctrl$time_h, ctrl$od, params)
    t0_res <- susi_detect_t0(trt$time_h, trt$od, params)
    ti_res <- if (is.na(t0_res$t0)) list(ti = NA_real_) else susi_detect_ti(trt$time_h, trt$od, t0_res$t0, params)
    susi_val <- susi_calc_susi(time_h, ctrl$od[match(time_h, ctrl$time_h)], trt$od[match(time_h, trt$time_h)],
                                t0_res$t0, ti_res$ti, tc_res$tc)
    return(susi_plot_curve(time_h, ctrl$od[match(time_h, ctrl$time_h)], trt$od[match(time_h, trt$time_h)],
                            t0_res$t0, ti_res$ti, tc_res$tc,
                            title = condition, subtitle = sprintf("overall mean | SusI = %.3f", susi_val)))
  }

  # biological_replicates
  bio_reps <- sort(unique(od_long$bio_rep[od_long$condition == condition]))
  if (!is.null(bio_rep)) bio_reps <- intersect(bio_reps, bio_rep)
  panels <- list()
  for (br in bio_reps) {
    ctrl <- wide_curve(control_all[control_all$bio_rep == br, ])
    trt  <- wide_curve(od_long[od_long$condition == condition & od_long$bio_rep == br, ])
    if (nrow(ctrl) == 0 || nrow(trt) == 0) next
    tc_res <- susi_detect_tc(ctrl$time_h, ctrl$od, params)
    t0_res <- susi_detect_t0(trt$time_h, trt$od, params)
    ti_res <- if (is.na(t0_res$t0)) list(ti = NA_real_) else susi_detect_ti(trt$time_h, trt$od, t0_res$t0, params)
    c_al <- ctrl$od[match(time_h, ctrl$time_h)]; p_al <- trt$od[match(time_h, trt$time_h)]
    susi_val <- susi_calc_susi(time_h, c_al, p_al, t0_res$t0, ti_res$ti, tc_res$tc)
    panels[[length(panels) + 1]] <- susi_plot_curve(time_h, c_al, p_al, t0_res$t0, ti_res$ti, tc_res$tc,
                                                      title = paste0(condition, " (bio_rep ", br, ")"),
                                                      subtitle = sprintf("SusI = %.3f", susi_val))
  }
  if (length(panels) == 0) stop("No matching control/treated curves found.", call. = FALSE)
  if (length(panels) == 1) return(panels[[1]])
  patchwork::wrap_plots(panels, ncol = min(3, length(panels)))
}

#' One-page diagnostic grid across every condition in the dataset
#'
#' Produces one panel per condition (pooled/overall-mean curves, for a
#' compact one-page overview; use [plot_condition()] for per-replicate
#' detail on any condition that looks off) and optionally saves it to a
#' file. Works for any number/labeling of conditions -- nothing here is
#' specific to MOI dilution series.
#'
#' @inheritParams run_susi
#' @param output_file Optional path (`.pdf` or `.png`) to save the grid to.
#' @param ncol Number of panel columns. Default: `min(3, n_conditions)`.
#' @return A combined `patchwork` object (also saved to `output_file` if given).
#' @export
plot_all_conditions <- function(file_path,
                                 od_sheet = "R", map_sheet = "name",
                                 well_col = "num", condition_col = "sample", bio_rep_col = "bio_rep",
                                 tech_reps_per_bio_rep = 12,
                                 control_label = NULL, conditions = NULL, exclude = NULL,
                                 time_limit_hours = NULL,
                                 params = susi_default_params(),
                                 output_file = NULL, ncol = NULL) {
  if (!is.null(time_limit_hours)) params$time_limit_hours <- time_limit_hours
  data <- read_plate_data(file_path, od_sheet, map_sheet, well_col, condition_col, bio_rep_col, tech_reps_per_bio_rep)
  cond_info <- susi_resolve_conditions(data$mapping[[condition_col]], control_label, conditions, exclude)
  control <- cond_info$control
  conds <- cond_info$conditions

  od_long <- data$od_long
  time_h <- data$time_h
  wide_curve <- function(df) stats::aggregate(od ~ time_h, data = df, FUN = mean, na.rm = TRUE)
  ctrl <- wide_curve(od_long[od_long$condition == control, ])
  tc_res <- susi_detect_tc(ctrl$time_h, ctrl$od, params)
  c_al <- ctrl$od[match(time_h, ctrl$time_h)]

  panels <- lapply(conds, function(cond) {
    trt <- wide_curve(od_long[od_long$condition == cond, ])
    t0_res <- susi_detect_t0(trt$time_h, trt$od, params)
    ti_res <- if (is.na(t0_res$t0)) list(ti = NA_real_) else susi_detect_ti(trt$time_h, trt$od, t0_res$t0, params)
    p_al <- trt$od[match(time_h, trt$time_h)]
    susi_val <- susi_calc_susi(time_h, c_al, p_al, t0_res$t0, ti_res$ti, tc_res$tc)
    susi_plot_curve(time_h, c_al, p_al, t0_res$t0, ti_res$ti, tc_res$tc,
                     title = cond, subtitle = sprintf("SusI = %.3f", susi_val))
  })

  if (is.null(ncol)) ncol <- min(3, length(panels))
  combined <- patchwork::wrap_plots(panels, ncol = ncol) +
    patchwork::plot_annotation(title = paste0("susiR diagnostic overview  |  control = '", control, "'"))

  if (!is.null(output_file)) {
    nr <- ceiling(length(panels) / ncol)
    ggplot2::ggsave(output_file, combined, width = 4 * ncol, height = 3.2 * nr, limitsize = FALSE)
    message("Diagnostic grid saved to ", output_file)
  }
  combined
}
