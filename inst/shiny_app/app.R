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
})

default_params <- susi_default_params()
`%||%` <- function(a, b) if (is.null(a) || identical(a, "")) b else a

# ============================================================================
# UI
# ============================================================================
ui <- fluidPage(
  titlePanel("susiR — phage lytic activity analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$details(
        tags$summary("Don't have an R / name sheet yet? Convert a raw plate-reader export"),
        br(),
        fileInput("raw_file", "Raw export from any plate reader (.xlsx)", accept = ".xlsx"),
        helpText("Auto-detects the data table (any sheet name, header row, well-ID spelling like A1/A:1, and time column/unit) and writes a converted file in the format below."),
        actionButton("convert_raw", "Convert", class = "btn-secondary", width = "100%"),
        verbatimTextOutput("raw_convert_summary"),
        uiOutput("raw_convert_download_ui")
      ),
      hr(),
      fileInput("file", "Plate-reader workbook (.xlsx)", accept = ".xlsx"),
      helpText("Requires two sheets: an OD time-course (time \u00d7 well) and a sample-mapping table (one row per well, with a condition label for each)."),

      checkboxInput("show_advanced", "Advanced: sheet / column names", value = FALSE),
      conditionalPanel(
        condition = "input.show_advanced",
        textInput("od_sheet", "OD sheet name", value = "R"),
        textInput("map_sheet", "Mapping sheet name", value = "name"),
        textInput("well_col", "Well-ID column", value = "num"),
        textInput("condition_col", "Condition-label column", value = "sample"),
        textInput("bio_rep_col", "Biological-replicate column (optional)", value = "bio_rep"),
        numericInput("tech_reps", "Fallback block size (wells per bio rep, if no bio_rep column)", value = 12, min = 1)
      ),

      hr(),
      uiOutput("control_ui"),
      uiOutput("exclude_ui"),
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
        tabPanel("Overview",
                 br(),
                 plotOutput("overview_plot", height = "700px"),
                 downloadButton("dl_overview", "Download overview (PNG)")),
        tabPanel("Combined view",
                 br(),
                 helpText("Every condition's mean curve on one plot, against the control."),
                 checkboxInput("combined_error_band", "Show \u00b1 SD across biological replicates", value = TRUE),
                 plotOutput("combined_plot", height = "550px"),
                 downloadButton("dl_combined", "Download combined view (PNG)")),
        tabPanel("Metrics summary",
                 br(),
                 helpText("SusI / VI / SupI / ti-tc across conditions: individual wells (circles), biological-replicate means (triangles), and the condition mean \u00b1 SD (black diamond and error bar), colour-coded by replicate."),
                 plotOutput("superplot", height = "650px"),
                 downloadButton("dl_superplot", "Download this view (PNG)")),
        tabPanel("Condition detail",
                 br(),
                 uiOutput("detail_condition_ui"),
                 plotOutput("detail_plot", height = "450px"),
                 downloadButton("dl_detail", "Download this plot (PNG)")),
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

  file_path <- reactive({
    req(input$file)
    input$file$datapath
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

  ## --- Main analysis, triggered only by the Run button -------------------
  result <- eventReactive(input$run, {
    withProgress(message = "Running susiR...", value = 0.3, {
      run_susi(file_path(),
               od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
               well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
               bio_rep_col = input$bio_rep_col %||% "bio_rep",
               tech_reps_per_bio_rep = input$tech_reps %||% 12,
               control_label = input$control_label, exclude = input$exclude_labels,
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
                           tech_reps_per_bio_rep = input$tech_reps %||% 12,
                           control_label = input$control_label, exclude = input$exclude_labels,
                           params = current_params())
    })
  })

  diag <- eventReactive(input$run, {
    diagnose_conditions(file_path(),
                         od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                         well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                         bio_rep_col = input$bio_rep_col %||% "bio_rep",
                         tech_reps_per_bio_rep = input$tech_reps %||% 12,
                         control_label = input$control_label, exclude = input$exclude_labels,
                         params = current_params(), verbose = FALSE)
  })

  output$overview_plot <- renderPlot({ req(overview()); overview() })

  combined <- eventReactive(input$run, {
    plot_combined_curves(file_path(),
                          od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                          well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                          bio_rep_col = input$bio_rep_col %||% "bio_rep",
                          tech_reps_per_bio_rep = input$tech_reps %||% 12,
                          control_label = input$control_label, exclude = input$exclude_labels,
                          time_limit_hours = current_params()$time_limit_hours,
                          show_error_band = isTRUE(input$combined_error_band))
  })
  output$combined_plot <- renderPlot({ req(combined()); combined() })

  superplot_obj <- eventReactive(input$run, {
    plot_metric_superplot(file_path(),
                           od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                           well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                           bio_rep_col = input$bio_rep_col %||% "bio_rep",
                           tech_reps_per_bio_rep = input$tech_reps %||% 12,
                           control_label = input$control_label, exclude = input$exclude_labels,
                           time_limit_hours = current_params()$time_limit_hours,
                           supi_window_hours = input$supi_window_hours,
                           params = current_params())
  })
  output$superplot <- renderPlot({ req(superplot_obj()); superplot_obj() })

  output$detail_condition_ui <- renderUI({
    req(result())
    selectInput("detail_condition", "Condition", choices = result()$conditions)
  })

  detail_plot_obj <- reactive({
    req(input$detail_condition)
    plot_condition(file_path(), condition = input$detail_condition,
                    od_sheet = input$od_sheet %||% "R", map_sheet = input$map_sheet %||% "name",
                    well_col = input$well_col %||% "num", condition_col = input$condition_col %||% "sample",
                    bio_rep_col = input$bio_rep_col %||% "bio_rep",
                    tech_reps_per_bio_rep = input$tech_reps %||% 12,
                    control_label = input$control_label,
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

  output$global_vi_text <- renderPrint({
    req(result())
    gv <- result()$global_vi
    if (isTRUE(gv$applicable)) {
      cat(sprintf("Global VI:   %.4f\nGlobal MV50: %s\n",
                   gv$VI, if (is.na(gv$MV50)) "not reached within tested MOI range" else sprintf("%.3g", gv$MV50)))
    } else {
      cat("Not applicable:\n", gv$note, "\n")
    }
  })

  output$dl_overview <- downloadHandler(
    filename = function() "susiR_overview.png",
    content = function(file) ggsave(file, overview(), width = 12, height = 8, dpi = 150)
  )
  output$dl_combined <- downloadHandler(
    filename = function() "susiR_combined_view.png",
    content = function(file) ggsave(file, combined(), width = 8, height = 5, dpi = 150)
  )
  output$dl_superplot <- downloadHandler(
    filename = function() "susiR_metrics_summary.png",
    content = function(file) ggsave(file, superplot_obj(), width = 11, height = 8, dpi = 150)
  )
  output$dl_detail <- downloadHandler(
    filename = function() paste0("susiR_", input$detail_condition, ".png"),
    content = function(file) ggsave(file, detail_plot_obj(), width = 9, height = 4, dpi = 150)
  )
  output$dl_excel <- downloadHandler(
    filename = function() "susiR_results.xlsx",
    content = function(file) {
      wb <- createWorkbook()
      addWorksheet(wb, "Summary"); writeData(wb, "Summary", result()$summary)
      addWorksheet(wb, "Diagnostics"); writeData(wb, "Diagnostics", diag())
      gv <- result()$global_vi
      addWorksheet(wb, "Global_VI")
      writeData(wb, "Global_VI", data.frame(VI = gv$VI, MV50 = gv$MV50, applicable = gv$applicable, note = gv$note %||% ""))
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
