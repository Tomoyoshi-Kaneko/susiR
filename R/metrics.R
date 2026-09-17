# ============================================================================
# metrics.R -- SusI, VI, SupI, ti/tc, computed from a control curve C(t), a
#               treated curve P(t), and detected time points t0, ti, tc.
# ============================================================================

#' Sustainability Index (SusI)
#'
#' SusI = \[ \eqn{\int_{t0}^{ti} C(t)dt - \int_{t0}^{ti} P(t)dt} \] /
#' \eqn{\int_{0}^{tc} C_{denom}(t)dt}, via trapezoidal integration.
#'
#' @param time_h Numeric vector of time points (hours), shared by `C` and `P`.
#' @param C Control (phage-free) OD curve used for the **numerator**
#'   (t0->ti window). For `calculation_method = "biological_replicates"`
#'   this is that specific biological replicate's own control curve.
#' @param P Phage-treated OD curve.
#' @param t0 Lysis-onset time (hours).
#' @param ti Resistance-emergence time (hours).
#' @param tc Control stationary-phase time (hours) (or, if already capped
#'   by a time limit, `effective_end_time`).
#' @param C_denom Control curve used for the **denominator**
#'   (0->tc window). Defaults to `C`. Pass the *pooled* (all biological
#'   replicates combined) control curve here when computing SusI per
#'   biological replicate, to match the legacy scripts: the denominator is
#'   one shared normalization constant per dataset, not one per replicate.
#' @return Numeric SusI value (unitless), or `NA` if inputs are insufficient.
#' @export
susi_calc_susi <- function(time_h, C, P, t0, ti, tc, C_denom = C) {
  if (any(is.na(c(t0, ti, tc)))) return(NA_real_)
  num <- susi_trapz(time_h, C, t0, ti) - susi_trapz(time_h, P, t0, ti)
  den <- susi_trapz(time_h, C_denom, 0, tc)
  if (is.na(den) || den == 0) return(NA_real_)
  num / den
}

#' Local Virulence Index contribution (VI), following Storms et al. (2020)
#'
#' Ported to match the legacy `calculate_vi_biological_replicates_v20()`:
#' integrates the *full* curve from `t = 0` up to `end_time` (the same
#' auto-detected control-stationary-phase time used elsewhere -- pass
#' `min(tc, time_limit_hours)`), **not** a short early-infection window.
#' An earlier version of this function used a fixed early window and
#' produced local-virulence values roughly half the size they should be;
#' this signature matches the validated legacy behavior.
#'
#' @inheritParams susi_calc_susi
#' @param end_time Upper bound of the integration window (hours) --
#'   normally `min(tc, time_limit_hours)`.
#' @return Numeric local-virulence value in \[0, 1\] (approximately), or `NA`.
#' @export
susi_calc_vi_local <- function(time_h, C, P, end_time) {
  auc_c <- susi_trapz(time_h, C, 0, end_time)
  auc_p <- susi_trapz(time_h, P, 0, end_time)
  if (is.na(auc_c) || auc_c == 0) return(NA_real_)
  (auc_c - auc_p) / auc_c
}

#' Global Virulence Index (VI) and MV50 across an MOI dilution series
#'
#' Fits the local-virulence-vs-MOI relationship (Storms et al., 2020) to
#' obtain a single VI summary value and the MOI (MV50) at which half-maximal
#' virulence is reached. Only applicable when conditions form a genuine MOI
#' dilution series -- check with [susi_check_moi_applicable()] first.
#'
#' @param moi_numeric Numeric MOI value per condition.
#' @param local_virulence Local virulence value per condition (same order).
#' @return A list: `VI`, `MV50` (both numeric, possibly `NA` if the fit fails).
#' @export
susi_calc_global_vi <- function(moi_numeric, local_virulence) {
  ok <- is.finite(moi_numeric) & is.finite(local_virulence) & moi_numeric > 0
  if (sum(ok) < 3) return(list(VI = NA_real_, MV50 = NA_real_))
  x <- log10(moi_numeric[ok]); y <- local_virulence[ok]
  o <- order(x); x <- x[o]; y <- y[o]

  # VI: area under the local-virulence-vs-log10(MOI) curve, normalized by
  # the maximum theoretically achievable area (rectangle of height 1 over
  # the observed log-MOI range) -- consistent with a bounded [0,1] index.
  auc <- susi_trapz(x, y, min(x), max(x))
  span <- max(x) - min(x)
  VI <- if (is.na(auc) || span <= 0) NA_real_ else auc / span

  # MV50: MOI at which local virulence first reaches half of its observed
  # maximum, by linear interpolation between bracketing points.
  half <- max(y, na.rm = TRUE) / 2
  MV50 <- NA_real_
  for (i in seq_len(length(x) - 1)) {
    if ((y[i] - half) * (y[i + 1] - half) <= 0 && y[i] != y[i + 1]) {
      frac <- (half - y[i]) / (y[i + 1] - y[i])
      MV50 <- 10^(x[i] + frac * (x[i + 1] - x[i]))
      break
    }
  }
  list(VI = VI, MV50 = MV50)
}

#' Suppression Index (SupI), following Kim et al. (2024)
#'
#' Percent reduction in area under the curve, relative to control, over a
#' fixed observation window (default 30 h, adjustable).
#'
#' @inheritParams susi_calc_susi
#' @param window_hours Length of the observation window (hours), from t = 0.
#' @return Numeric SupI value in percent, or `NA`.
#' @export
susi_calc_supi <- function(time_h, C, P, window_hours = 30) {
  auc_c <- susi_trapz(time_h, C, 0, window_hours)
  auc_p <- susi_trapz(time_h, P, 0, window_hours)
  if (is.na(auc_c) || auc_c == 0) return(NA_real_)
  (auc_c - auc_p) / auc_c * 100
}

#' Time-ratio metric (ti/tc)
#'
#' @param ti Resistance-emergence time (hours).
#' @param tc Control stationary-phase time (hours).
#' @return Numeric ratio, or `NA`.
#' @export
susi_calc_time_ratio <- function(ti, tc) {
  if (any(is.na(c(ti, tc))) || tc == 0) return(NA_real_)
  ti / tc
}
