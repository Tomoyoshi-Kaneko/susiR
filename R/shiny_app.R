# ============================================================================
# shiny_app.R -- launcher for the bundled Shiny GUI (inst/shiny_app/app.R)
# ============================================================================

#' Launch the susiR Shiny application
#'
#' Opens a browser-based interface implementing the full susiR workflow --
#' file upload, automatic condition/replicate detection, interactive
#' detection-parameter tuning with live diagnostic plots, and Excel/plot
#' download -- without requiring any R code to be written.
#'
#' @param ... Passed on to [shiny::runApp()] (e.g. `launch.browser = TRUE`,
#'   `port = 1234`).
#' @export
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The 'shiny' package is required to run the app. Install it with install.packages('shiny').", call. = FALSE)
  }
  app_dir <- system.file("shiny_app", package = "susiR")
  if (!nzchar(app_dir)) {
    # Development mode (susiR not yet installed as a package): assume the
    # current working directory is the susiR package root, exactly as used
    # throughout this project so far.
    app_dir <- file.path(getwd(), "inst", "shiny_app")
    if (!dir.exists(app_dir)) {
      stop("Could not locate the Shiny app directory ('inst/shiny_app'). ",
           "Please run this with your working directory set to the susiR package root.", call. = FALSE)
    }
  }
  shiny::runApp(app_dir, ...)
}
