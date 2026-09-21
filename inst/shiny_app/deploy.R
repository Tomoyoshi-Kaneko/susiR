# ============================================================================
# deploy.R -- deployment notes for the susiR Shiny app.
#
# ACTUAL WORKFLOW IN USE: GitHub -> Connect Cloud ("Publish from GitHub").
# This repo is connected to Posit Connect Cloud directly; every push to the
# connected branch (main) triggers Connect Cloud to rebuild and redeploy
# the app automatically -- usually within a minute or two. Most of the
# time, deploying an update is simply:
#
#   1. Edit files locally (R/ and/or inst/shiny_app/).
#   2. If the set of R packages the app needs changed (a new library() call
#      was added anywhere in R/ or inst/shiny_app/app.R), regenerate
#      manifest.json (see below) -- this is the one step that's easy to
#      forget, and the one thing Connect Cloud needs that a plain git push
#      doesn't automatically figure out on its own.
#   3. Push everything to GitHub (drag-and-drop upload is fine).
#   4. Wait for Connect Cloud to auto-redeploy, or trigger it manually from
#      the app's page on connect.posit.cloud if it doesn't pick it up.
#
# No local R session, browser-based login, or rsconnect deployApp() call is
# needed for this path -- it's driven entirely by what's in the repo.
# ============================================================================

## --- Regenerating manifest.json (only needed when dependencies change) ----
## Run this from R, from the susiR package root (same working directory
## used throughout the README), with every package the app uses actually
## installed locally first:
##   install.packages(c("shiny", "DT", "ggplot2", "patchwork", "openxlsx",
##                       "dplyr", "tidyr", "xml2", "readxl", "rhandsontable",
##                       "svglite", "rsconnect"))
rsconnect::writeManifest(
  appDir   = "inst/shiny_app",
  appFiles = c(list.files("inst/shiny_app", full.names = FALSE),
               file.path("..", "..", "R", list.files("R")))
)
## This writes inst/shiny_app/manifest.json, listing every package (and the
## installed version of each) the app needs. Commit and push this file
## alongside your other changes -- Connect Cloud reads it to know what to
## install when it rebuilds the app. If you forget this step after adding
## a new library() call, the rebuilt app will fail with a "package not
## found" style error rather than silently ignoring the new code.

## ----------------------------------------------------------------------
## ALTERNATIVE (not currently used): deploying directly from R via
## rsconnect, without going through GitHub at all. This is documented here
## in case the GitHub path is ever unavailable, but as of this writing
## rsconnect::connectCloudUser() has a known bug in at least the current
## release that makes it fail immediately after a successful browser login
## (a "serverName must be a single string, not NULL" error, reproducible
## even in fresh sessions and even after upgrading to the GitHub dev build
## of rsconnect). The GitHub-based workflow above sidesteps this entirely
## and is what this project actually uses.
##
## If you want to try it again later (e.g. after an rsconnect bug fix):
##   install.packages("rsconnect")
##   rsconnect::connectCloudUser()   # opens a browser to sign in
##   rsconnect::accounts()           # confirm it registered
##   rsconnect::deployApp(
##     appDir      = "inst/shiny_app",
##     appFiles    = c(list.files("inst/shiny_app", full.names = FALSE),
##                     file.path("..", "..", "R", list.files("R"))),
##     appName     = "susiR",
##     appTitle    = "susiR -- phage lytic activity analysis",
##     forceUpdate = TRUE
##   )
## ----------------------------------------------------------------------
