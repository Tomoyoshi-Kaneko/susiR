# ============================================================================
# app.R -- susiR Shiny application
#
# Full workflow without writing any R code: upload a plate-reader workbook,
# susiR auto-detects the experimental structure (conditions, control,
# biological replicates), lets you tune detection parameters with a live
# diagnostic plot, then computes SusI/VI/SupI/ti-tc for every condition and
# lets you download the results (Excel) and figures (PNG).
# ============================================================================

## --- Load susiR core functions --------------------------------------------
## When this app is launched via susiR::run_app() from an *installed*
## package, the core functions are already attached via library(susiR)
## below. In development mode (this project, pre-installation), shiny sets
## the working directory to this app's own folder while it runs, so the
## package's R/ folder is reachable two levels up.
.core_dir <- file.path("..", "..", "R")
if (dir.exists(.core_dir)) {
  for (.f in list.files(.core_dir, pattern = "\\.R$", full.names = TRUE)) source(.f)
} else if (requireNamespace("susiR", quietly = TRUE)) {
  library(susiR)
}

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(xml2)
  if (requireNamespace("readxl", quietly = TRUE)) library(readxl)
  if (requireNamespace("rhandsontable", quietly = TRUE)) library(rhandsontable)
  if (requireNamespace("svglite", quietly = TRUE)) library(svglite)
})

default_params <- susi_default_params()
`%||%` <- function(a, b) if (is.null(a) || identical(a, "")) b else a

# ============================================================================
# UI
# ============================================================================
ui <- fluidPage(
  titlePanel(
    div(
      style = "display: flex; align-items: center;",
      tags$img(src = "logo.png", height = "60px", style = "margin-right: 15px;"),
      div(
        h2("susiR", style = "margin: 0; font-weight: bold; color: #8E1728;"), # 早稲田エンジ色
        h5("phage lytic activity analysis", style = "margin: 0; color: #555;")
      )
    ),
    windowTitle = "susiR — phage lytic activity analysis"
  ),
  tags$head(tags$style(HTML("
    summary.susi-toggle {
      cursor: pointer; padding: 8px 12px; margin-bottom: 8px;
      background-color: #e8eef7; border: 1px solid #9fb3d1; border-radius: 5px;
      font-weight: bold; display: block;
    }
    summary.susi-toggle:hover { background-color: #d8e2f2; }
    summary.susi-toggle::after { content: '  \u25be click to expand / collapse'; font-weight: normal; color: #555; }
    details[open] > summary.susi-toggle::after { content: '  \u25b4 click to collapse'; }
  "))),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$details(
        tags$summary(class = "susi-toggle", "Don't have an R / name sheet yet? Convert a raw plate-reader export"),
        br(),
        fileInput("raw_file", "Raw export from any plate reader (.xlsx)", accept = ".xlsx"),
        helpText("Auto-detects the data table (any sheet name, header row, well-ID spelling like A1/A:1, and time column/unit) and writes a converted file in the format below."),
        actionButton("convert_raw", "Convert", class = "btn-secondary", width = "100%"),
        verbatimTextOutput("raw_convert_summary"),
        uiOutput("raw_convert_download_ui")
      ),
      tags$details(
        tags$summary(class = "susi-toggle", "\U0001F5BC\uFE0F Export settings: image size & PNG/SVG format"),
        br(),
        numericInput("export_width", "Width (inches)", value = 10, min = 1, step = 0.5),
        numericInput("export_height", "Height (inches)", value = 8, min = 1, step = 0.5),
        helpText("Increase height especially for Overview/Metrics summary when there are many panels and curves look compressed."),
        selectInput("export_format", "Image format", choices = c("PNG" = "png", "SVG (vector, editable)" = "svg"))
      ),
      hr(),
      fileInput("file", "Plate-reader workbook (.xlsx)", accept = ".xlsx"),
      helpText("Requires two sheets: an OD time-course (time \u00d7 well) and a sample-mapping table (one row per well, with a condition label for each)."),

      checkboxInput("show_advanced", "Advanced: sheet / column names", value = FALSE),
      conditionalPanel(
        condition = "input.show_advanced",
        uiOutput("advanced_settings_ui")
      ),

      checkboxInput("use_host", "Analyze in separate groups (e.g. different hosts/strains sharing this sheet)", value = FALSE),
      conditionalPanel("input.use_host",
        uiOutput("host_col_ui"),
        helpText("Each group's control is auto-detected by matching common control-like labels (\"ct\", \"control\", \"no phage\", \"blank\", etc., or any label containing \"host\") -- this won't catch every naming convention (e.g. \"cntrl\", or a plain strain name used as the control with no such word in it). Override any group below if needed."),
        uiOutput("control_override_ui"),
        helpText("The Overview tab always shows every group at once; the selector below only picks which group the other plot tabs (Combined view, Metrics summary, Condition detail) focus on."),
        uiOutput("host_select_ui"),
        conditionalPanel("input.run > 0", strong("Control used per group (after Run):"), verbatimTextOutput("group_controls_info"))
      ),

      hr(),
      conditionalPanel("!input.use_host", uiOutput("control_ui")),
      uiOutput("exclude_ui"),
      textInput("condition_order_input", "Display order (optional, comma-separated)",
                value = "", placeholder = "e.g. 10^0, 10^-1, 10^-2"),
      textInput("condition_colors_input", "Custom colors (optional, label=color, comma-separated)",
                value = "", placeholder = "e.g. 10^0=red, 10^-1=steelblue"),
      helpText("Anything not named in \"Display order\" keeps its usual order, placed after the ones you did name. Anything not named in \"Custom colors\" keeps the default palette. For colors, common names work directly (red, blue, green, orange, purple, steelblue, ...); for anything more specific, a hex code like #1b9e77 works too."),
      selectInput("calculation_method", "Aggregation level",
                  choices = c("biological_replicates", "individual_wells", "overall_mean"),
                  selected = "biological_replicates"),
      checkboxInput("use_time_limit", "Limit analysis to a maximum time (h)", value = TRUE),
      conditionalPanel("input.use_time_limit", numericInput("time_limit_hours", NULL, value = 24, min = 1)),
      numericInput("supi_window_hours", "SupI window (h)", value = 30, min = 1),

      hr(),
      h4("Detection parameters"),
      sliderInput("smooth_window", "Smoothing window (points)", min = 1, max = 11, value = default_params$smooth_window, step = 2),
      sliderInput("t0_min_time", "Earliest allowed t0 (h)", min = 0, max = 5, value = default_params$t0_min_time, step = 0.05),
      sliderInput("t0_sustainability", "t0 search window length (h)", min = 1, max = 30, value = default_params$t0_sustainability, step = 1),
      sliderInput("n_forward_points_t0", "t0: forward-looking points", min = 2, max = 30, value = default_params$n_forward_points_t0, step = 1),
      sliderInput("threshold_percentage_t0", "t0: confirmation fraction", min = 0.1, max = 1, value = default_params$threshold_percentage_t0, step = 0.05),
      sliderInput("min_lysis_time", "Min. time from t0 to look for ti (h)", min = 0, max = 10, value = default_params$min_lysis_time, step = 0.25),
      sliderInput("n_forward_points_ti", "ti: forward-looking points", min = 2, max = 30, value = default_params$n_forward_points_ti, step = 1),
      sliderInput("threshold_percentage_ti", "ti: confirmation fraction", min = 0.1, max = 1, value = default_params$threshold_percentage_ti, step = 0.05),
      sliderInput("min_increase_threshold", "ti: min. OD increase", min = 0, max = 0.1, value = default_params$min_increase_threshold, step = 0.001),
      sliderInput("tc_min_time", "Earliest allowed tc (h)", min = 0, max = 10, value = default_params$tc_min_time, step = 0.25),
      sliderInput("tc_slope_threshold", "tc: slope threshold", min = 0, max = 0.2, value = default_params$tc_slope_threshold, step = 0.005),

      hr(),
      actionButton("run", "Run analysis", class = "btn-primary", width = "100%")
    ),

    mainPanel(
      width = 8,
      tabsetPanel(
        id = "tabs",
        tabPanel("Build sample map (visual)",
                 br(),
                 helpText("An alternative to preparing a 'name' sheet in Excel: click-and-drag across the grid below to select a range of wells, type a label, and press Ctrl+Enter to fill the whole selection at once (standard spreadsheet behavior). Wells greyed out aren't present in the OD sheet currently selected under Advanced settings. Only wells you fill in the 'Sample / condition' grid are included -- leave a well blank to exclude it."),
                 textOutput("plate_sheet_info"),
                 tabsetPanel(
                   tabPanel("Sample / condition (required)", br(), rHandsontableOutput("grid_sample")),
                   tabPanel("Group (optional -- e.g. host/strain)", br(), rHandsontableOutput("grid_group")),
                   tabPanel("Biological replicate (optional)", br(), rHandsontableOutput("grid_biorep"))
                 ),
                 br(),
                 actionButton("build_mapping", "Build mapping from grid and use it below", class = "btn-primary"),
                 verbatimTextOutput("build_mapping_summary")),
        tabPanel("Overview",
                 br(),
                 plotOutput("overview_plot", height = "700px"),
                 downloadButton("dl_overview", "Download overview")),
        tabPanel("Combined view",
                 br(),
                 helpText("Every condition's mean curve on one plot, against the control."),
                 checkboxInput("combined_error_band", "Show \u00b1 SD across biological replicates", value = TRUE),
                 plotOutput("combined_plot", height = "550px"),
                 downloadButton("dl_combined", "Download combined view")),
        tabPanel("Metrics summary",
                 br(),
                 helpText("SusI / VI / SupI / ti-tc across conditions: individual wells (circles), biological-replicate means (triangles), and the condition mean \u00b1 SD (black diamond and error bar), colour-coded by replicate."),
                 plotOutput("superplot", height = "650px"),
                 downloadButton("dl_superplot", "Download this view"),
                 hr(),
                 helpText("Simple bar-chart view of the same numbers (respects \"Display order\" and \"Custom colors\" from the sidebar):"),
                 plotOutput("bar_chart", height = "400px"),
                 downloadButton("dl_bar_chart", "Download this view")),
        tabPanel("Condition detail",
                 br(),
                 uiOutput("detail_condition_ui"),
                 plotOutput("detail_plot", height = "450px"),
                 downloadButton("dl_detail", "Download this plot")),
        tabPanel("Results table",
                 br(),
                 DTOutput("summary_table"),
                 br(),
                 downloadButton("dl_excel", "Download results (Excel)")),
        tabPanel("Diagnostics",
                 br(),
                 helpText("Per-well/per-replicate detection outcome. Anything other than 'success' is worth a look on the Condition detail tab."),
                 DTOutput("diag_table")),
        tabPanel("Global VI / MV50",
                 br(),
                 verbatimTextOutput("global_vi_text"))
      )
    )
  )
)

# ============================================================================
# Server
# ============================================================================
server <- function(input, output, session) {

  raw_converted_path <- reactiveVal(NULL)

  observeEvent(input$convert_raw, {
    req(input$raw_file)
    out_path <- file.path(tempdir(), paste0("converted_", tools::file_path_sans_ext(basename(input$raw_file$name)), ".xlsx"))
    result <- tryCatch({
      withCallingHandlers(
        susi_import_raw(input$raw_file$datapath, output_path = out_path),
        message = function(m) {
          output$raw_convert_summary <- renderPrint(cat(conditionMessage(m)))
          invokeRestart("muffleMessage")
        }
      )
    }, error = function(e) {
      output$raw_convert_summary <- renderPrint(cat("Could not auto-convert this file:\n", conditionMessage(e)))
      NULL
    })
    if (!is.null(result)) raw_converted_path(out_path)
  })

  output$raw_convert_download_ui <- renderUI({
    req(raw_converted_path())
    downloadButton("dl_raw_converted", "Download converted file", class = "btn-secondary", width = "100%")
  })
  output$dl_raw_converted <- downloadHandler(
    filename = function() "converted_for_susiR.xlsx",
    content = function(file) file.copy(raw_converted_path(), file)
  )

  visual_mapping_path <- reactiveVal(NULL)

  ## Uploading a new file always starts fresh: any mapping built visually
  ## for a previous file (and its grid contents) no longer applies.
  observeEvent(input$file, {
    visual_mapping_path(NULL)
    grid_sample_data(.well_grid_template())
    grid_group_data(.well_grid_template())
    grid_biorep_data(.well_grid_template())
  }, ignoreInit = TRUE)

  file_path <- reactive({
    if (!is.null(visual_mapping_path())) return(visual_mapping_path())
    req(input$file)
    input$file$datapath
  })

  ## --- Sheet-name and host-column discovery (dropdowns, not free text) --
  sheet_names <- reactive({
    req(file_path())
    tryCatch(.sheet_names_any(file_path()), error = function(e) character(0))
  })

  ## Sheet/column settings, consolidated into one dynamic block so they can
  ## be locked together: once a mapping has been built visually, every
  ## field here is fixed to match it (num/sample/bio_rep, sheets R/name) --
  ## editing any of them would silently stop matching what was built,
  ## which is exactly the confusion this avoids. Re-upload the file (which
  ## resets the visual mapping, above) to edit these again.
  output$advanced_settings_ui <- renderUI({
    if (!is.null(visual_mapping_path())) {
      return(tagList(
        helpText(strong("Fixed by the visual mapping editor:"), "OD sheet 'R', mapping sheet 'name', well column 'num', condition column 'sample', bio_rep column 'bio_rep'. Re-upload the file above to edit these directly instead."),
        selectInput("od_sheet", NULL, choices = "R", selected = "R"),
        selectInput("map_sheet", NULL, choices = "name", selected = "name")
      ))
    }
    sn <- sheet_names()
    od_sel <- if ("R" %in% sn) "R" else if (length(sn) > 0) sn[1] else NULL
    map_sel <- if ("name" %in% sn) "name" else if (length(sn) > 1) sn[2] else if (length(sn) > 0) sn[1] else NULL
    tagList(
      selectInput("od_sheet", "OD sheet name", choices = sn, selected = od_sel),
      selectInput("map_sheet", "Mapping sheet name", choices = sn, selected = map_sel),
      textInput("well_col", "Well-ID column", value = "num"),
      textInput("condition_col", "Condition-label column", value = "sample"),
      textInput("bio_rep_col", "Biological-replicate column (optional)", value = "bio_rep"),
      numericInput("tech_reps", "Max wells per biological replicate (fallback)", value = 12, min = 1),
      helpText("Only used when there's no bio_rep column above. If a condition has this many wells or fewer, they're all automatically treated as ONE biological replicate -- leave this at the default if you're unsure or don't have separate biological replicates (e.g. 4 wells per condition, all technical replicates: the default of 12 already gives the correct answer, 1 replicate of 4, with nothing to change).")
    )
  })

  mapping_columns <- reactive({
    req(file_path(), input$map_sheet)
    tryCatch(names(openxlsx::read.xlsx(file_path(), sheet = input$map_sheet)), error = function(e) character(0))
  })
  output$host_col_ui <- renderUI({
    cols <- mapping_columns()
    selectInput("host_col_choice", "Grouping column (e.g. host, strain -- whatever splits this sheet into separate analyses, each with its own control)", choices = cols, selected = if ("host" %in% cols) "host" else cols[1])
  })

  host_col_effective <- reactive({ if (isTRUE(input$use_host)) input$host_col_choice else NULL })

  host_values <- reactive({
    req(input$use_host, input$host_col_choice, file_path(), input$map_sheet)
    tryCatch({
      map <- openxlsx::read.xlsx(file_path(), sheet = input$map_sheet)
      v <- unique(as.character(map[[input$host_col_choice]]))
      v[!is.na(v) & nzchar(trimws(v))]
    }, error = function(e) character(0))
  })
  output$host_select_ui <- renderUI({
    hv <- host_values()
    selectInput("selected_host", "Group to display (plots other than Overview)", choices = hv, selected = if (length(hv) > 0) hv[1] else NULL)
  })

  ## Per-group control override: lets a user pick, group by group, which of
  ## that group's own labels is its control -- for naming conventions the
  ## automatic detection doesn't recognize (anything besides "ct",
  ## "control", "no phage", "blank", or a label containing "host").
  group_condition_labels <- reactive({
    hv <- host_values()
    req(length(hv) > 0, file_path(), input$map_sheet)
    cc <- input$condition_col %||% "sample"
    tryCatch({
      map <- openxlsx::read.xlsx(file_path(), sheet = input$map_sheet)
      stats::setNames(lapply(hv, function(h) sort(unique(as.character(map[[cc]][map[[input$host_col_choice]] == h])))), hv)
    }, error = function(e) list())
  })

  output$control_override_ui <- renderUI({
    gl <- group_condition_labels()
    if (length(gl) == 0) return(NULL)
    tagList(
      lapply(seq_along(gl), function(i) {
        selectInput(paste0("ctrl_override_", i), names(gl)[i],
                    choices = c("(auto-detect)", gl[[i]]), selected = "(auto-detect)")
      })
    )
  })

  group_control_overrides <- reactive({
    gl <- group_condition_labels()
    if (length(gl) == 0) return(NULL)
    vals <- vapply(seq_along(gl), function(i) input[[paste0("ctrl_override_", i)]] %||% "(auto-detect)", character(1))
    names(vals) <- names(gl)
    vals <- vals[vals != "(auto-detect)"]
    if (length(vals) == 0) return(NULL)
    vals
  })

  ## --- Visual 96-well plate mapping editor --------------------------------
  ## The sheet to read OD data from for the plate editor: prefer the user's
  ## explicit selection, but fall back to a sensible default (first sheet
  ## found, or "R" if present) rather than silently doing nothing -- an
  ## empty/unset dropdown value should never block reading the file.
  plate_od_sheet <- reactive({
    if (!is.null(input$od_sheet) && nzchar(input$od_sheet)) return(input$od_sheet)
    sn <- sheet_names()
    if ("R" %in% sn) "R" else if (length(sn) > 0) sn[1] else NULL
  })

  plate_valid_wells <- reactive({
    req(file_path(), plate_od_sheet())
    tryCatch({
      od <- openxlsx::read.xlsx(file_path(), sheet = plate_od_sheet())
      setdiff(names(od), "Time")
    }, error = function(e) character(0))
  })

  output$plate_sheet_info <- renderText({
    if (!is.null(visual_mapping_path())) {
      return("Mapping already built from this grid. Editing the grid below and clicking the button again will rebuild and replace it; re-upload the file above to start over completely.")
    }
    paste0("Determining which wells exist by reading OD sheet: '", plate_od_sheet() %||% "(upload a file first)",
           "' (this only checks which wells have data, to grey out the rest below -- it does not read any sample names).")
  })

  .well_grid_template <- function() {
    m <- matrix("", nrow = 8, ncol = 12)
    rownames(m) <- LETTERS[1:8]
    colnames(m) <- as.character(1:12)
    as.data.frame(m, stringsAsFactors = FALSE, check.names = FALSE)
  }

  grid_sample_data <- reactiveVal(.well_grid_template())
  grid_group_data  <- reactiveVal(.well_grid_template())
  grid_biorep_data <- reactiveVal(.well_grid_template())

  .make_plate_hot <- function(df, valid_wells) {
    rh <- rhandsontable::rhandsontable(df, rowHeaderWidth = 50, width = 760, height = 320) |>
      rhandsontable::hot_cols(colWidths = 55)
    for (r in 1:8) {
      for (cc in 1:12) {
        wid <- paste0(LETTERS[r], cc)
        if (!(wid %in% valid_wells)) rh <- rhandsontable::hot_cell(rh, r, cc, readOnly = TRUE)
      }
    }
    rh
  }

  output$grid_sample <- rhandsontable::renderRHandsontable(.make_plate_hot(grid_sample_data(), plate_valid_wells()))
  output$grid_group  <- rhandsontable::renderRHandsontable(.make_plate_hot(grid_group_data(),  plate_valid_wells()))
  output$grid_biorep <- rhandsontable::renderRHandsontable(.make_plate_hot(grid_biorep_data(), plate_valid_wells()))

  observeEvent(input$grid_sample, { grid_sample_data(rhandsontable::hot_to_r(input$grid_sample)) })
  observeEvent(input$grid_group,  { grid_group_data(rhandsontable::hot_to_r(input$grid_group)) })
  observeEvent(input$grid_biorep, { grid_biorep_data(rhandsontable::hot_to_r(input$grid_biorep)) })

  observeEvent(input$build_mapping, {
    req(file_path(), plate_od_sheet())
    vw <- plate_valid_wells()
    s <- grid_sample_data(); g <- grid_group_data(); b <- grid_biorep_data()
    rows <- list()
    for (r in 1:8) {
      for (cc in 1:12) {
        wid <- paste0(LETTERS[r], cc)
        if (!(wid %in% vw)) next
        samp <- s[r, cc]
        if (is.na(samp) || !nzchar(trimws(samp))) next
        grp <- g[r, cc]; br <- b[r, cc]
        rows[[length(rows) + 1]] <- data.frame(
          num = wid, sample = trimws(samp),
          host = if (!is.na(grp) && nzchar(trimws(grp))) trimws(grp) else NA_character_,
          bio_rep = if (!is.na(br) && nzchar(trimws(br))) trimws(br) else NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows) == 0) {
      output$build_mapping_summary <- renderPrint(cat("No wells filled in yet in the 'Sample / condition' grid -- nothing to build."))
      return(invisible())
    }
    map_df <- dplyr::bind_rows(rows)
    has_host <- !all(is.na(map_df$host))
    has_biorep <- !all(is.na(map_df$bio_rep))
    if (!has_host) map_df$host <- NULL
    if (!has_biorep) map_df$bio_rep <- NULL

    od <- openxlsx::read.xlsx(file_path(), sheet = plate_od_sheet())
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "R"); openxlsx::writeData(wb, "R", od)
    openxlsx::addWorksheet(wb, "name"); openxlsx::writeData(wb, "name", map_df)
    out_path <- file.path(tempdir(), "susiR_visual_mapping.xlsx")
    openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
    visual_mapping_path(out_path)

    output$build_mapping_summary <- renderPrint({
      cat(sprintf("Built mapping for %d well(s). This is now the active file below -- OD sheet 'R' and mapping sheet 'name' are selected automatically.", nrow(map_df)))
      if (has_host) cat("\nGroup column written as 'host' -- check \"Analyze in separate groups\" above and pick 'host' if you want to use it.")
      if (has_biorep) cat("\nBiological-replicate column 'bio_rep' detected and will be used automatically.")
    })
  })

  current_params <- reactive({
    p <- default_params
    for (nm in setdiff(names(p), "time_limit_hours")) {
      if (!is.null(input[[nm]])) p[[nm]] <- input[[nm]]
    }
    p$time_limit_hours <- if (isTRUE(input$use_time_limit)) input$time_limit_hours else NULL
    p
  })

  ## --- Discover condition labels as soon as a file is uploaded (and
  ## whenever the advanced sheet/column settings change), so the
  ## control/exclude selectors can be populated without a full run. --------
  mapping_labels <- reactive({
    req(file_path())
    tryCatch({
      map <- openxlsx::read.xlsx(file_path(), sheet = input$map_sheet %||% "name")
      col <- input$condition_col %||% "sample"
      validate(need(col %in% names(map), paste0("Column '", col, "' not found in mapping sheet.")))
      unique(as.character(map[[col]]))
    }, error = function(e) {
      validate(paste("Could not read mapping sheet:", conditionMessage(e)))
    })
  })

  guessed_control <- reactive({
    labs <- mapping_labels()
    tryCatch(susi_resolve_conditions(labs, control_label = NULL)$control, error = function(e) labs[1])
  })

  output$control_ui <- renderUI({
    labs <- mapping_labels()
    selectInput("control_label", "Control condition", choices = labs, selected = guessed_control())
  })

  output$exclude_ui <- renderUI({
    labs <- mapping_labels()
    selectInput("exclude_labels", "Exclude condition(s) (optional)", choices = labs, selected = NULL, multiple = TRUE)
  })

  condition_order_parsed <- reactive({
    txt <- trimws(input$condition_order_input %||% "")
    if (!nzchar(txt)) return(NULL)
    parts <- trimws(strsplit(txt, ",")[[1]])
    parts[nzchar(parts)]
  })

  condition_colors_parsed <- reactive({
    txt <- trimws(input$condition_colors_input %||% "")
    if (!nzchar(txt)) return(NULL)
    pairs <- strsplit(txt, ",")[[1]]
    kv <- lapply(pairs, function(p) trimws(strsplit(p, "=")[[1]]))
    kv <- kv[vapply(kv, length, integer(1)) == 2]
    if (length(kv) == 0) return(NULL)
    stats::setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))
  })

  ## --- Main analysis, triggered only by the Run button -------------------
  result <- eventReactive(input$run, {
    withProgress(message = "Running susiR...", value = 0.3, {
      run_susi(file_path(),
               od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
               well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
               bio_rep_col = input$bio_rep_col %||% "bio_rep",
               host_col = host_col_effective(),
               tech_reps_per_bio_rep = input$tech_reps %||% 12,
               control_label = if (isTRUE(input$use_host)) group_control_overrides() else input$control_label,
               exclude = input$exclude_labels, condition_order = condition_order_parsed(),
               calculation_method = input$calculation_method,
               supi_window_hours = input$supi_window_hours,
               params = current_params(), verbose = FALSE)
    })
  })

  overview <- eventReactive(input$run, {
    withProgress(message = "Building overview plot...", value = 0.7, {
      plot_all_conditions(file_path(),
                           od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                           well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                           bio_rep_col = input$bio_rep_col %||% "bio_rep",
                           host_col = host_col_effective(),
                           tech_reps_per_bio_rep = input$tech_reps %||% 12,
                           control_label = if (isTRUE(input$use_host)) group_control_overrides() else input$control_label,
                           exclude = input$exclude_labels, condition_order = condition_order_parsed(),
                           params = current_params())
    })
  })

  diag <- eventReactive(input$run, {
    diagnose_conditions(file_path(),
                         od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                         well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                         bio_rep_col = input$bio_rep_col %||% "bio_rep",
                         host_col = host_col_effective(),
                         tech_reps_per_bio_rep = input$tech_reps %||% 12,
                         control_label = if (isTRUE(input$use_host)) group_control_overrides() else input$control_label,
                         exclude = input$exclude_labels, condition_order = condition_order_parsed(),
                         params = current_params(), verbose = FALSE)
  })

  output$overview_plot <- renderPlot({ req(overview()); overview() })

  combined <- eventReactive(input$run, {
    plot_combined_curves(file_path(),
                          od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                          well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                          bio_rep_col = input$bio_rep_col %||% "bio_rep",
                          host_col = host_col_effective(), host = if (isTRUE(input$use_host)) input$selected_host else NULL,
                          tech_reps_per_bio_rep = input$tech_reps %||% 12,
                          control_label = if (isTRUE(input$use_host)) group_control_overrides() else input$control_label,
                          exclude = input$exclude_labels, condition_order = condition_order_parsed(),
                          time_limit_hours = current_params()$time_limit_hours,
                          show_error_band = isTRUE(input$combined_error_band),
                          condition_colors = condition_colors_parsed())
  })
  output$combined_plot <- renderPlot({ req(combined()); combined() })

  superplot_obj <- eventReactive(input$run, {
    plot_metric_superplot(file_path(),
                           od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                           well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                           bio_rep_col = input$bio_rep_col %||% "bio_rep",
                           host_col = host_col_effective(), host = if (isTRUE(input$use_host)) input$selected_host else NULL,
                           tech_reps_per_bio_rep = input$tech_reps %||% 12,
                           control_label = if (isTRUE(input$use_host)) group_control_overrides() else input$control_label,
                           exclude = input$exclude_labels, condition_order = condition_order_parsed(),
                           time_limit_hours = current_params()$time_limit_hours,
                           supi_window_hours = input$supi_window_hours,
                           params = current_params())
  })
  output$superplot <- renderPlot({ req(superplot_obj()); superplot_obj() })

  bar_chart_obj <- eventReactive(input$run, {
    if (isTRUE(input$use_host)) {
      req(input$selected_host)
      r <- result()
      sub_summary <- r$summary[r$summary$host == input$selected_host, ]
      sub_result <- list(summary = sub_summary, conditions = r$conditions[[input$selected_host]])
      plot_metrics_summary(sub_result, condition_colors = condition_colors_parsed())
    } else {
      plot_metrics_summary(result(), condition_colors = condition_colors_parsed())
    }
  })
  output$bar_chart <- renderPlot({ req(bar_chart_obj()); bar_chart_obj() })

  output$detail_condition_ui <- renderUI({
    req(result())
    choices <- if (isTRUE(input$use_host)) result()$conditions[[input$selected_host]] else result()$conditions
    selectInput("detail_condition", "Condition", choices = choices)
  })

  detail_plot_obj <- reactive({
    req(input$detail_condition)
    plot_condition(file_path(), condition = input$detail_condition,
                    od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                    well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                    bio_rep_col = input$bio_rep_col %||% "bio_rep",
                    host_col = host_col_effective(), host = if (isTRUE(input$use_host)) input$selected_host else NULL,
                    tech_reps_per_bio_rep = input$tech_reps %||% 12,
                    control_label = if (isTRUE(input$use_host)) group_control_overrides() else input$control_label,
                    calculation_method = if (input$calculation_method == "individual_wells") "biological_replicates" else input$calculation_method,
                    params = current_params())
  })
  output$detail_plot <- renderPlot({ detail_plot_obj() })

  output$summary_table <- renderDT({
    req(result())
    datatable(result()$summary, options = list(pageLength = 15, scrollX = TRUE)) |>
      formatRound(intersect(c("SusI","VI_local","SupI","ti_tc","t0","ti","tc","SusI_se","VI_local_se","SupI_se","ti_tc_se"), names(result()$summary)), 3)
  })

  output$diag_table <- renderDT({
    req(diag())
    datatable(diag(), options = list(pageLength = 20, scrollX = TRUE)) |>
      formatStyle("status", backgroundColor = styleEqual(
        c("success", "t0_detection_failed", "ti_not_detected"),
        c("#e8f5e9", "#ffebee", "#fff8e1")
      ))
  })

  output$group_controls_info <- renderPrint({
    req(result())
    ctrl <- result()$control
    if (isTRUE(input$use_host) && !is.null(names(ctrl))) {
      for (h in names(ctrl)) cat(h, "->", ctrl[[h]], "\n")
    } else {
      cat(ctrl, "\n")
    }
  })

  output$global_vi_text <- renderPrint({
    req(result())
    gv <- result()$global_vi
    if (isTRUE(input$use_host)) {
      for (h in names(gv)) {
        cat("Host:", h, "\n")
        if (isTRUE(gv[[h]]$applicable)) {
          cat(sprintf("  Global VI:   %.4f\n  Global MV50: %s\n",
                       gv[[h]]$VI, if (is.na(gv[[h]]$MV50)) "not reached within tested MOI range" else sprintf("%.3g", gv[[h]]$MV50)))
        } else {
          cat("  Not applicable:", gv[[h]]$note, "\n")
        }
        cat("\n")
      }
    } else if (isTRUE(gv$applicable)) {
      cat(sprintf("Global VI:   %.4f\nGlobal MV50: %s\n",
                   gv$VI, if (is.na(gv$MV50)) "not reached within tested MOI range" else sprintf("%.3g", gv$MV50)))
    } else {
      cat("Not applicable:\n", gv$note, "\n")
    }
  })

  .export_ext <- function() input$export_format %||% "png"
  .save_export <- function(file, plot) {
    fmt <- .export_ext()
    w <- input$export_width %||% 10
    h <- input$export_height %||% 8
    if (fmt == "svg") {
      ggsave(file, plot, width = w, height = h, device = svglite::svglite)
    } else {
      ggsave(file, plot, width = w, height = h, dpi = 150)
    }
  }

  output$dl_overview <- downloadHandler(
    filename = function() paste0("susiR_overview.", .export_ext()),
    content = function(file) .save_export(file, overview())
  )
  output$dl_combined <- downloadHandler(
    filename = function() paste0("susiR_combined_view.", .export_ext()),
    content = function(file) .save_export(file, combined())
  )
  output$dl_superplot <- downloadHandler(
    filename = function() paste0("susiR_metrics_summary.", .export_ext()),
    content = function(file) .save_export(file, superplot_obj())
  )
  output$dl_bar_chart <- downloadHandler(
    filename = function() paste0("susiR_bar_chart.", .export_ext()),
    content = function(file) .save_export(file, bar_chart_obj())
  )
  output$dl_detail <- downloadHandler(
    filename = function() paste0("susiR_", input$detail_condition, ".", .export_ext()),
    content = function(file) .save_export(file, detail_plot_obj())
  )
  output$dl_excel <- downloadHandler(
    filename = function() "susiR_results.xlsx",
    content = function(file) {
      wb <- createWorkbook()
      addWorksheet(wb, "Summary"); writeData(wb, "Summary", result()$summary)
      addWorksheet(wb, "Diagnostics"); writeData(wb, "Diagnostics", diag())
      gv <- result()$global_vi
      gv_df <- if (isTRUE(input$use_host)) {
        dplyr::bind_rows(lapply(names(gv), function(h) data.frame(
          host = h, VI = gv[[h]]$VI, MV50 = gv[[h]]$MV50, applicable = gv[[h]]$applicable, note = gv[[h]]$note %||% ""
        )))
      } else {
        data.frame(VI = gv$VI, MV50 = gv$MV50, applicable = gv$applicable, note = gv$note %||% "")
      }
      addWorksheet(wb, "Global_VI")
      writeData(wb, "Global_VI", gv_df)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
