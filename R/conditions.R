# ============================================================================
# conditions.R -- generalized condition/control resolution. No fixed set of
#                 MOI labels is ever assumed: conditions can be MOI values,
#                 phage names, or any other label the user has used.
# ============================================================================

.default_control_patterns <- c("ct", "control", "no phage", "no-phage", "nophage",
                                "phage-free", "phage free", "uninfected", "blank")

#' Identify the control condition and the set of conditions to analyze
#'
#' No assumption is made about how many conditions exist or what they are
#' called -- MOI dilutions, different phages at a single fixed MOI, phage
#' cocktails, anything. Every distinct label in the mapping table's
#' condition column is treated as one condition; exactly one of them is
#' identified as the phage-free control, and the rest become the conditions
#' analyzed by [run_susi()].
#'
#' @param condition_labels Character vector of condition labels as they
#'   appear (repeated per well) in the mapping table.
#' @param control_label Optional explicit control label. If given, it must
#'   be one of `condition_labels`.
#' @param conditions Optional explicit vector of condition labels to
#'   analyze (also fixes their order, e.g. for plotting). Defaults to every
#'   non-control label, in order of first appearance.
#' @param exclude Optional vector of labels to drop entirely (e.g. a
#'   `"free"` / phage-only-no-bacteria well type that is neither the control
#'   nor a condition of interest).
#' @return A list with `control` (character scalar) and `conditions`
#'   (character vector, analysis order).
#' @export
susi_resolve_conditions <- function(condition_labels, control_label = NULL,
                                     conditions = NULL, exclude = NULL) {
  labels <- unique(condition_labels)

  if (is.null(control_label)) {
    hit <- labels[tolower(trimws(labels)) %in% .default_control_patterns]
    if (length(hit) == 1) {
      control_label <- hit
      message("Control condition auto-detected as '", control_label, "'.")
    } else if (length(hit) == 0) {
      stop("Could not auto-detect the control condition among: ",
           paste(labels, collapse = ", "),
           ". Please pass `control_label` explicitly.", call. = FALSE)
    } else {
      stop("Multiple labels look like a control condition (",
           paste(hit, collapse = ", "), "). Please pass `control_label` explicitly.", call. = FALSE)
    }
  } else if (!control_label %in% labels) {
    stop("`control_label` = '", control_label, "' does not match any condition label in the data (",
         paste(labels, collapse = ", "), ").", call. = FALSE)
  }

  remaining <- setdiff(labels, control_label)
  if (!is.null(exclude)) remaining <- setdiff(remaining, exclude)

  if (is.null(conditions)) {
    conditions <- remaining[order(match(remaining, unique(condition_labels)))]
  } else {
    unknown <- setdiff(conditions, labels)
    if (length(unknown) > 0) {
      stop("`conditions` includes label(s) not present in the data: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }

  list(control = control_label, conditions = conditions)
}

#' Decide whether a set of condition labels supports a dose-response
#' (MOI-based) global VI/MV50 fit
#'
#' Global VI (per Storms et al. 2020) is only meaningful when conditions
#' represent a dilution series of a single phage against a single host, so
#' that "local virulence" can be plotted/fit against MOI. If conditions are
#' e.g. different phages tested at one fixed MOI, this does not apply, and
#' [run_susi()] skips the global fit (while still reporting local, per-
#' condition VI, SusI, SupI and ti/tc, which never require this assumption).
#'
#' @param conditions Character vector of condition labels.
#' @return A list with `applicable` (logical) and, if applicable,
#'   `moi_numeric` (named numeric vector, condition label -> MOI).
#' @export
susi_check_moi_applicable <- function(conditions) {
  moi <- susi_parse_moi_numeric(conditions)
  if (all(!is.na(moi)) && length(unique(moi)) == length(moi) && length(moi) >= 3) {
    list(applicable = TRUE, moi_numeric = moi)
  } else {
    list(applicable = FALSE, moi_numeric = moi)
  }
}
