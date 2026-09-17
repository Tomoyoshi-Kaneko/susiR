# ============================================================================
# utils.R -- small numerical / parsing helpers shared across susiR
# ============================================================================

#' Centered moving average
#'
#' @param x Numeric vector.
#' @param window_size Window width (odd integers behave most symmetrically).
#' @return Smoothed numeric vector, same length as `x`.
#' @keywords internal
susi_moving_average <- function(x, window_size = 3) {
  n <- length(x)
  if (n < window_size || window_size <= 1) return(x)
  half <- floor(window_size / 2)
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }
  out
}

#' Trapezoidal-rule definite integral of a curve, restricted to [from, to]
#'
#' Linearly interpolates the curve at the `from`/`to` boundaries if they do
#' not fall exactly on an observed time point, so integration windows do not
#' need to be snapped to the sampling grid.
#'
#' @param time Numeric vector of time points (strictly increasing).
#' @param value Numeric vector of curve values, same length as `time`.
#' @param from Lower integration bound (defaults to `min(time)`).
#' @param to Upper integration bound (defaults to `max(time)`).
#' @return A single numeric value (area under the curve), or `NA` if the
#'   window is degenerate or data are insufficient.
#' @keywords internal
susi_trapz <- function(time, value, from = NULL, to = NULL) {
  ok <- !is.na(time) & !is.na(value)
  time <- time[ok]; value <- value[ok]
  if (length(time) < 2) return(NA_real_)
  o <- order(time); time <- time[o]; value <- value[o]

  if (is.null(from)) from <- time[1]
  if (is.null(to))   to   <- time[length(time)]
  if (!is.finite(from) || !is.finite(to) || to <= from) return(NA_real_)
  from <- max(from, time[1])
  to   <- min(to, time[length(time)])
  if (to <= from) return(NA_real_)

  interp <- function(t0) {
    if (t0 <= time[1]) return(value[1])
    if (t0 >= time[length(time)]) return(value[length(time)])
    stats::approx(time, value, xout = t0)$y
  }

  keep <- time > from & time < to
  tt <- c(from, time[keep], to)
  vv <- c(interp(from), value[keep], interp(to))

  sum(diff(tt) * (utils::head(vv, -1) + utils::tail(vv, -1)) / 2)
}

#' Sort key for plate well IDs such as "A1", "A12", "H9"
#'
#' Orders first by row letter(s), then numerically by column number, so
#' `"A2" < "A12" < "B1"` (a plain string sort would put `"A12"` before
#' `"A2"`).
#'
#' @param wells Character vector of well IDs.
#' @return Integer rank vector suitable for `order()`.
#' @keywords internal
susi_well_sort_key <- function(wells) {
  m <- regmatches(wells, regexec("^([A-Za-z]+)([0-9]+)$", wells))
  letters_part <- vapply(m, function(x) if (length(x) == 3) x[2] else NA_character_, character(1))
  nums_part    <- vapply(m, function(x) if (length(x) == 3) as.integer(x[3]) else NA_integer_, integer(1))
  # Convert letter prefix to a base-26-ish rank (handles A..Z, AA.. if ever needed)
  letter_rank <- vapply(letters_part, function(s) {
    if (is.na(s)) return(NA_integer_)
    v <- utf8ToInt(toupper(s)) - utf8ToInt("A") + 1
    as.integer(Reduce(function(acc, d) acc * 26 + d, v, accumulate = FALSE))
  }, integer(1))
  order(letter_rank, nums_part, na.last = TRUE)
}

#' Does a character vector look like MOI-style numeric labels?
#'
#' Recognizes plain numbers ("1", "0.5") and `"10^-5"`-style exponent
#' notation, either of which is common for MOI dilution series.
#'
#' @param labels Character vector of condition labels.
#' @return Named numeric vector (NA for labels that do not parse), same
#'   length/order as `labels`.
#' @keywords internal
susi_parse_moi_numeric <- function(labels) {
  parse_one <- function(s) {
    s <- trimws(s)
    if (grepl("^10\\^-?[0-9]+(\\.[0-9]+)?$", s)) {
      exp_txt <- sub("^10\\^", "", s)
      return(10^as.numeric(exp_txt))
    }
    if (grepl("^-?[0-9]+(\\.[0-9]+)?([eE]-?[0-9]+)?$", s)) {
      return(as.numeric(s))
    }
    NA_real_
  }
  out <- vapply(labels, parse_one, numeric(1))
  names(out) <- labels
  out
}
