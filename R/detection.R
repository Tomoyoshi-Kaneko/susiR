# ============================================================================
# detection.R -- automatic detection of t0 (lysis onset), ti (resistance
#                 emergence), and tc (control stationary phase). Operates on
#                 a single (time, OD) curve; the same code runs regardless
#                 of how many/which conditions are present in the dataset.
#
#                 Ported to match the validated legacy algorithm from the
#                 original analysis scripts (detect_lysis_start_time_peak_based(),
#                 detect_resistance_emergence_forward(), determine_tc()) after
#                 a numeric-validation pass found the from-scratch
#                 reimplementation differed from it in several concrete ways
#                 (bounded t0 search window, true local-maximum requirement,
#                 monotonic-fallback time, sustained tc confirmation, and the
#                 ti search-start / confirmation rule) -- see README.
# ============================================================================

#' Default detection parameters
#'
#' Defaults follow the values used for the bundled T1/MG1655 example dataset
#' (see the accompanying Application Note draft, Table S1, and the legacy
#' `process_phage_data_v20()` quick-start call). All are exposed so they can
#' be retuned per-dataset without touching any code.
#'
#' @return A named list of detection parameters.
#' @export
susi_default_params <- function() {
  list(
    smooth_window            = 3,
    t0_min_time               = 0.25,
    t0_sustainability          = 20,
    n_forward_points_t0       = 10,
    threshold_percentage_t0   = 0.6,
    n_forward_points_ti       = 7,
    threshold_percentage_ti   = 0.6,
    min_increase_threshold    = 0.008,
    min_lysis_time            = 2.0,
    tc_slope_threshold         = 0.03,
    tc_min_time                = 1.0,
    time_limit_hours           = NULL
  )
}

#' Detect lysis onset (t0) on a single OD curve
#'
#' Ported to match the validated legacy algorithm
#' (`detect_lysis_start_time_peak_based()` in the original analysis
#' scripts): searches for a genuine **local maximum** -- OD at or above
#' both neighbors, strictly above at least one -- within a *bounded* window
#' `[t0_min_time, t0_min_time + t0_sustainability]`, confirmed by a
#' sustained decline over the next `n_forward_points_t0` observations (at
#' least `threshold_percentage_t0` of them lower). If the curve is already
#' declining from the very first eligible point (>= 80% of the first up to
#' 10 pairwise comparisons are decreasing), `t0` is set to the *first
#' observed time point overall* (not `t0_min_time`) and flagged
#' `"monotonic"`, matching the legacy behavior. If no peak is found within
#' the bounded window, `t0 = NA` is returned (the legacy scripts' internal
#' fallback value of 0 is never actually used downstream by their own
#' calling code -- a failed well is excluded from calculation there too, so
#' `NA` is the faithful equivalent at this function's boundary).
#'
#' @param time_h Numeric vector of time points (hours).
#' @param od Numeric vector of OD values, same length as `time_h`.
#' @param params List of detection parameters, see [susi_default_params()].
#' @return A list: `t0` (numeric, or `NA` on failure), `method` (`"peak"`,
#'   `"monotonic"`, or `"failed"`).
#' @export
susi_detect_t0 <- function(time_h, od, params = susi_default_params()) {
  o <- order(time_h); time_h <- time_h[o]; od <- od[o]
  od_s <- susi_moving_average(od, params$smooth_window)
  n <- length(time_h)

  start_idx <- which(time_h >= params$t0_min_time)[1]
  if (is.na(start_idx)) return(list(t0 = NA_real_, method = "failed"))

  ## --- monotonic-decrease check (legacy: fixed first up-to-10 points) ----
  if (start_idx >= 1) {
    check_length <- min(10, n - start_idx)
    if (check_length >= 1) {
      decreases <- 0
      for (i in start_idx:(start_idx + check_length - 1)) {
        if (!is.na(od_s[i]) && !is.na(od_s[i + 1]) && od_s[i + 1] < od_s[i]) decreases <- decreases + 1
      }
      if (decreases / check_length >= 0.8) {
        return(list(t0 = time_h[1], method = "monotonic"))
      }
    }
  }

  ## --- bounded local-maximum search --------------------------------------
  search_end_time <- params$t0_min_time + params$t0_sustainability
  search_end_idx <- utils::tail(which(time_h <= search_end_time), 1)
  if (length(search_end_idx) == 0 || start_idx >= search_end_idx) return(list(t0 = NA_real_, method = "failed"))

  last_i <- search_end_idx - params$n_forward_points_t0
  if (last_i < start_idx) return(list(t0 = NA_real_, method = "failed"))

  for (i in start_idx:last_i) {
    if (i < 2 || i >= n) next
    cur <- od_s[i]; prev <- od_s[i - 1]; nxt <- od_s[i + 1]
    if (any(is.na(c(cur, prev, nxt)))) next

    is_local_max <- (cur >= prev) && (cur >= nxt) && ((cur > prev) || (cur > nxt))
    if (!is_local_max) next

    hi <- min(i + params$n_forward_points_t0, n)
    fwd <- od_s[(i + 1):hi]
    fwd <- fwd[!is.na(fwd)]
    if (length(fwd) < ceiling(params$n_forward_points_t0 * 0.5)) next

    if (mean(fwd < cur) >= params$threshold_percentage_t0) {
      return(list(t0 = time_h[i], method = "peak"))
    }
  }
  list(t0 = NA_real_, method = "failed")
}

#' Detect resistance emergence (ti) on a single OD curve
#'
#' Ported to match the validated legacy algorithm
#' (`detect_resistance_emergence_forward()`): the search for sustained
#' regrowth starts at whichever is *later* of (a) the time of the actual
#' minimum OD observed after `t0`, or (b) `t0 + min_lysis_time` -- not
#' simply `t0 + min_lysis_time` as such. A candidate point is confirmed
#' when at least `threshold_percentage_ti` of the next `n_forward_points_ti`
#' observations are higher than it *and* their **average** increase is at
#' least `min_increase_threshold` (not a per-point count).
#'
#' A successfully detected `ti` is capped at `params$time_limit_hours` (if
#' set); a *failed* detection falls back to `fallback_limit` (the caller
#' should pass `min(tc, time_limit_hours)` here, matching the legacy
#' scripts' `effective_end_time` -- not `time_limit_hours` alone).
#'
#' @inheritParams susi_detect_t0
#' @param t0 Lysis-onset time (hours), typically from [susi_detect_t0()].
#' @param fallback_limit Value to return when no regrowth is detected.
#'   Defaults to `params$time_limit_hours` (or the last observed time
#'   point) if not supplied; pass `min(tc, time_limit_hours)` to match the
#'   legacy behavior exactly.
#' @return A list: `ti` (numeric), `detected` (logical).
#' @export
susi_detect_ti <- function(time_h, od, t0, params = susi_default_params(), fallback_limit = NULL) {
  o <- order(time_h); time_h <- time_h[o]; od <- od[o]
  od_s <- susi_moving_average(od, params$smooth_window)
  n <- length(time_h)

  cap <- if (!is.null(params$time_limit_hours)) min(max(time_h), params$time_limit_hours) else max(time_h)
  limit <- if (!is.null(fallback_limit)) fallback_limit else cap

  after_t0 <- which(time_h > t0 & !is.na(od_s))
  if (length(after_t0) == 0) return(list(ti = limit, detected = FALSE))
  min_idx <- after_t0[which.min(od_s[after_t0])]
  min_od_time <- time_h[min_idx]

  search_start_time <- max(min_od_time, t0 + params$min_lysis_time)
  start_idx <- which(time_h >= search_start_time)[1]
  end_idx <- n - params$n_forward_points_ti
  if (is.na(start_idx) || start_idx > end_idx) return(list(ti = limit, detected = FALSE))

  for (i in start_idx:end_idx) {
    cur <- od_s[i]
    if (is.na(cur)) next
    fwd <- od_s[(i + 1):(i + params$n_forward_points_ti)]
    fwd <- fwd[!is.na(fwd)]
    if (length(fwd) < ceiling(params$n_forward_points_ti * 0.5)) next

    pct_higher <- mean(fwd > cur)
    if (pct_higher >= params$threshold_percentage_ti) {
      avg_increase <- mean(fwd - cur)
      if (!is.na(avg_increase) && avg_increase >= params$min_increase_threshold) {
        ti_time <- time_h[i]
        return(list(ti = min(ti_time, cap), detected = TRUE))
      }
    }
  }
  list(ti = limit, detected = FALSE)
}

#' Detect control stationary phase (tc)
#'
#' Ported to match the validated legacy algorithm (`determine_tc()`):
#' requires a **sustained** low growth rate, not a single point below
#' `tc_slope_threshold` -- once a candidate point's slope drops below the
#' threshold, at least `max(3, ceiling(0.6 * n_check))` of the next
#' `n_check` (up to 5) points must *also* be below threshold before `tc` is
#' accepted. This avoids one noisy dip triggering an early `tc`.
#'
#' @inheritParams susi_detect_t0
#' @return A list: `tc` (numeric), `detected` (logical).
#' @export
susi_detect_tc <- function(time_h, od, params = susi_default_params()) {
  o <- order(time_h); time_h <- time_h[o]; od <- od[o]
  n <- length(time_h)
  if (n < 2) return(list(tc = max(time_h), detected = FALSE))

  slopes <- numeric(n - 1)
  for (i in seq_len(n - 1)) {
    dt <- time_h[i + 1] - time_h[i]
    dOD <- od[i + 1] - od[i]
    slopes[i] <- if (!is.na(dt) && !is.na(dOD) && dt > 0) dOD / dt else NA_real_
  }
  below <- slopes <= params$tc_slope_threshold & !is.na(slopes)

  for (i in which(below)) {
    if (time_h[i] < params$tc_min_time) next
    n_check <- min(5, length(slopes) - i + 1)
    if (n_check < 3) next
    subsequent <- i:(i + n_check - 1)
    n_low <- sum(below[subsequent], na.rm = TRUE)
    required <- max(3, ceiling(n_check * 0.6))
    if (n_low >= required) return(list(tc = time_h[i], detected = TRUE))
  }
  list(tc = max(time_h), detected = FALSE)
}
