# susiR <img src="man/figures/logo.png" align="right" height="150" />

**susiR** is an open-source R package and interactive Shiny application for standardizing, automating, and code-free quantification of bacteriophage lytic dynamics from microplate growth-curve data.

susiR computes four complementary indices of bacteriophage lytic activity from optical density (OD) time-course data:

- **SusI** (Sustainability Index) — duration and depth of lysis suppression, from lysis onset (`t0`) to resistance emergence (`ti`)
- **VI** (Virulence Index) — early infection dynamics, following Storms et al. (2020)
- **SupI** (Suppression Index) — overall growth suppression over a defined observation window, following Kim et al. (2024)
- **ti/tc** — a time-ratio metric relating suppression duration to the control culture's normal growth cycle

susiR is the companion software for the Sustainability Index manuscript (Kaneko et al., in preparation). Unlike original custom analysis scripts, susiR makes no assumptions about experimental design: the number of conditions, their labels, and the number and layout of biological replicates are inferred from a user-supplied sample-mapping table at run time, allowing the same pipeline to handle diverse plate layouts without modifying any code.

---

## Try susiR first — no programming required

If you are not familiar with R or programming, **you do not need to install anything to try susiR**.

### New to susiR? Start with the visual guide

A short slide-deck walkthrough (with screenshots for every step) is available at [`docs/HowToUse.pdf`](docs/HowToUse.pdf).

### Use the web application

Open the susiR web application in your browser:

** [Launch susiR](https://01a0aedc-8749-7813-9276-04c501bdd47d.share.connect.posit.cloud/)**

The web application provides a point-and-click workflow for uploading data, checking the plate layout, adjusting analysis parameters, viewing diagnostic plots, and downloading results.

### Start with the example data

The easiest way to understand the workflow is to start with the bundled **T1 example dataset**:

[`inst/extdata/T1_rawdata.xlsx`](inst/extdata/T1_rawdata.xlsx)

You can use the example data directly in the web application where the sample/example-data options are provided.

The T1 workbook is also a useful **template for preparing your own input file**. If you have data from another experiment, you can create a workbook that follows the same basic structure as `T1_rawdata.xlsx`, or use the application's plate-map editor and raw-format importer where appropriate.

> **In short:** if you just want to try the software, open the web app and start with the T1 example. If you want to analyze your own experiment, you can either prepare an Excel file in the same format as T1 or use the tools provided in the app to help prepare it.

---

## What can susiR do?

### Key features

- **Flexible replicate count.** Biological replicates are not fixed at 3. An explicit `bio_rep` column in the mapping sheet is used exactly as given — 1, 2, 4, or any other number, arranged however the wells were laid out.

- **Flexible replicate layout.** Because replicate membership can be declared per well, replicates do not need to occupy a fixed, contiguous block of plate columns. Without a `bio_rep` column, susiR falls back to a positional convention (consecutive blocks of 12 wells per condition, in plate order).

- **Flexible condition labels.** Conditions are not assumed to be an MOI dilution series. A sheet can contain several phages tested at one fixed MOI, or any other labeling scheme. Every distinct label in the mapping sheet is treated as a condition. The MOI-dependent global VI/MV50 dose-response fit is computed only when condition labels can be interpreted as a numeric MOI series; local, per-condition SusI/VI/SupI/ti-tc are still computed for other labeling schemes.

- **Import from different plate readers.** Raw exports with different sheet names, header rows, well-ID formats, and time units can be converted into the format expected by susiR.

- **Visual plate-map editor.** The Shiny application provides an interactive alternative to preparing the `name` sheet manually.

- **Grouped analysis.** Multiple hosts or strains can be analyzed separately when they are present in the same workbook and each group has its own control.

- **Point-and-click interface.** The bundled Shiny application provides the full workflow — upload, condition/replicate detection, parameter tuning, diagnostic plots, and results download — without writing R code.

---

## Input data

susiR expects an Excel workbook with two main sheets:

1. an **OD data sheet**, normally named `R`
2. a **sample-mapping sheet**, normally named `name`

The exact sheet and column names can be configured when using the R functions.

### OD sheet

The default OD sheet is named `R`.

| Column | Meaning |
|---|---|
| `Time` | Measurement time |
| Other columns | One column for each well |

`Time` can be either:

- a plain number interpreted as elapsed seconds, or
- a genuine Excel time/duration cell, such as a cell formatted as `[h]:mm`.

`read_plate_data()` automatically inspects the workbook's cell formatting when `time_unit = "auto"` (the default), so Excel duration values are not accidentally interpreted as seconds.

### Mapping sheet

The default mapping sheet is named `name`. It contains one row per well.

| Column | Required? | Meaning |
|---|---|---|
| `well` | Yes | Well ID; must match a column in the OD sheet |
| `sample` | Yes | Condition label; any text is allowed |
| `bio_rep` | No | Explicit biological-replicate ID |

If `bio_rep` is omitted, susiR uses the positional fallback: consecutive blocks of 12 wells per condition, in plate order.

If `bio_rep` is present, it should be filled in for every well belonging to the relevant condition.

Column names can be changed with `well_col`, `condition_col`, and `bio_rep_col`.

### Preparing your own file

For users who are unfamiliar with programming, the simplest approach is to use the bundled T1 workbook as a template:

[`inst/extdata/T1_rawdata.xlsx`](inst/extdata/T1_rawdata.xlsx)

You can copy its structure and replace the OD measurements and sample labels with those from your own experiment.

You do **not** have to reproduce the original T1 experimental design. The important point is to preserve the relationship between the OD well columns and the corresponding rows in the mapping sheet.

---

## Importing data from other plate readers

Raw exports from plate-reader software often do not match the input format above. Different instruments may use different sheet names, header rows, well-ID formats, and time units.

`susi_import_raw()` can convert many such exports into a susiR-compatible workbook:

```r
susi_import_raw(
  "my_raw_export.xlsx",
  output_path = "converted.xlsx"
)
```

For example, the importer can recognize well IDs such as:

- `A1`
- `A:1`
- `A01`
- well IDs embedded in longer labels such as `"Sample0001 (B02)"`

It searches workbook sheets for a plausible time row and well columns, removes unrelated metadata columns, and writes a workbook containing:

- an `R` sheet with the converted OD/time data
- a `name` sheet containing the detected well IDs, ready for sample assignment

The importer recognizes common English and Japanese time labels such as `time` / `時間`, `sec` / `秒`, and `min` / `分`.

The printed import summary should always be checked before using the converted file. If the match is too uncertain, the importer fails explicitly rather than silently guessing.

The importer has been validated against exports from:

- BioTek/Agilent Epoch2
- Promega
- Thermo Fisher SkanIt, including Japanese-language exports

If an export is too unusual for automatic detection, you can specify the appropriate sheet with `sheet =`, or convert the file manually to match the bundled T1 example.

The Shiny application also provides an **"Import a raw plate-reader export"** workflow.

---

## Using susiR from R

The web application is the recommended starting point for users who do not use R.

For R users, susiR can be used directly from the repository without installing it as a package:

```r
for (f in list.files("R", full.names = TRUE)) {
  source(f)
}
```

A standard package workflow using `devtools::load_all()` or `R CMD build` is also possible.

---

## Basic usage

### T1 example

The bundled T1 dataset does not contain an explicit `bio_rep` column, so susiR uses the positional fallback.

```r
res <- run_susi(
  "inst/extdata/T1_rawdata.xlsx",
  control_label = "Ct",
  exclude = "free",
  calculation_method = "biological_replicates",
  time_limit_hours = 24
)

res$summary
res$global_vi
```

`res$summary` contains SusI, local VI, SupI, and ti/tc results by condition, including biological-replicate summaries where applicable.

`res$global_vi` contains the global VI/MV50 analysis when the condition labels form an applicable numeric MOI series.

### Custom data with explicit biological replicates

If your experiment has an irregular replicate structure, add a `bio_rep` column to the mapping sheet. No changes to the analysis code are required.

```r
res2 <- run_susi(
  "my_plate.xlsx",
  control_label = "Ct"
)
```

### Several phages at one fixed MOI

Conditions do not have to be an MOI dilution series.

```r
res3 <- run_susi(
  "cocktail_screen.xlsx",
  control_label = "Ct"
)

res3$global_vi$note
```

In this situation, the global VI/MV50 dose-response analysis is skipped because a numeric MOI series is not available, while the local SusI/VI/SupI/ti-tc calculations remain available for each condition.

### Diagnose conditions before running the full analysis

For fast parameter tuning and quality control:

```r
diag <- diagnose_conditions(
  "inst/extdata/T1_rawdata.xlsx",
  control_label = "Ct",
  exclude = "free",
  output_excel = "T1_diagnosis.xlsx"
)
```

---

## Plotting

Functions in `R/plotting.R` and `R/plotting_summary.R` generate diagnostic and summary figures. These functions require `ggplot2` and `patchwork`.

### Per-condition diagnostic plot

```r
plot_condition(
  "inst/extdata/T1_rawdata.xlsx",
  condition = "10^0",
  control_label = "Ct",
  time_limit_hours = 24
)
```

### All conditions

```r
plot_all_conditions(
  "inst/extdata/T1_rawdata.xlsx",
  control_label = "Ct",
  exclude = "free",
  time_limit_hours = 24,
  output_file = "diagnostic_overview.png"
)
```

### Combined growth curves

```r
plot_combined_curves(
  "inst/extdata/T1_rawdata.xlsx",
  control_label = "Ct",
  exclude = "free",
  time_limit_hours = 24
)
```

The plot shows each condition's mean growth curve, with mean ± SD across biological replicates.

### Index summary

```r
res <- run_susi(
  "inst/extdata/T1_rawdata.xlsx",
  control_label = "Ct",
  exclude = "free",
  calculation_method = "biological_replicates",
  time_limit_hours = 24
)

plot_metrics_summary(res)
```

### Well- and replicate-level superplot

```r
plot_metric_superplot(
  "inst/extdata/T1_rawdata.xlsx",
  control_label = "Ct",
  exclude = "free",
  time_limit_hours = 24
)

plot_metric_superplot(
  "inst/extdata/T1_rawdata.xlsx",
  metrics = "SusI",
  control_label = "Ct",
  exclude = "free",
  time_limit_hours = 24
)
```

`plot_metric_superplot()` displays individual wells, biological-replicate means, and condition-level means together, making within- and between-replicate variability visible.

The plotting functions adapt automatically to the condition labels. A numeric MOI series is displayed using the corresponding MOI ordering; other labeling schemes fall back to categorical plots.

`plot_metric_superplot()` runs `run_susi()` internally at both the `individual_wells` and `biological_replicates` levels, so it may take longer than the other plotting functions.

---

## Shiny application

The repository contains a point-and-click Shiny application at:

```text
inst/shiny_app/app.R
```

It covers the main workflow:

1. Upload an Excel file
2. Review the detected conditions and controls
3. Adjust detection parameters if necessary
4. Run the analysis
5. Inspect diagnostic and summary plots
6. Download result tables and figures

The application also includes:

- raw plate-reader import
- visual plate-map editing
- grouped analysis
- interactive parameter tuning
- Excel and PNG/SVG export

### Run locally

The Shiny application requires R.

Install the additional packages:

```r
install.packages(c("shiny", "DT", "patchwork"))
```

Then either:

- open `inst/shiny_app/app.R` in RStudio and click **Run App**, or
- from the package root, run:

```r
source("R/shiny_app.R")
run_app()
```

### Use the hosted application

If you do not want to install R locally, use the hosted application:

**[Launch susiR Web Application](https://01a0aedc-8749-7813-9276-04c501bdd47d.share.connect.posit.cloud/)**

---

## Deployment

The Shiny application can also be deployed to Posit Connect Cloud.

See:

```text
inst/shiny_app/deploy.R
```

for deployment instructions using either the `rsconnect` R package or publication from a GitHub repository.

---

## Validation

SusI and VI have been checked against the exact values reported in the companion manuscript for the T1/*E. coli* MG1655 dataset at MOI = 1:

| Index | susiR | Manuscript |
|---|---:|---:|
| SusI | 0.664 | 0.664 |
| VI (local, MOI = 1) | 0.86 | 0.84 |
| VI (global) | 0.84 | 0.84 |

The local VI difference reflects the current implementation and should not be interpreted as exact numerical identity with the manuscript's reported local value.

ti/tc and SupI were verified by re-deriving their calculation formulas directly from the original analysis scripts rather than against a single reported value; both are implemented to match those formulas.

The remaining five MOI conditions in the dataset have not been checked against exact tabulated manuscript values. The manuscript reports precise values only at MOI = 1, with the full dilution series shown graphically rather than tabulated. A comparison against that figure showed close agreement across the T1 series.

---

## Implementation notes

A few calculation details are worth stating explicitly for reproducibility.

### Control stationary-phase time (`tc`)

`tc` is computed once per dataset, from all control wells pooled together, and reused throughout:

- as the denominator bound for SusI
- as the integration bound for VI
- as the fallback endpoint for `ti`

SusI's numerator (the `t0` → `ti` window) uses each biological replicate's own control curve; only the denominator uses the pooled control curve.

### Virulence Index (VI)

VI is integrated from `t = 0` to the same control stationary-phase time (`tc`) used for SusI, following Storms et al. (2020), rather than over a fixed window to the end of the observation period.

### Suppression Index (SupI)

The SupI integration window defaults to `time_limit_hours` rather than a fixed 30 hours, matching how the original pipeline was run for the companion manuscript.

### Reported variability

`SusI_se`, `VI_local_se`, `SupI_se`, and `ti_tc_se` in `run_susi()$summary`, together with the error bars in `plot_metric_superplot()`, represent the standard error across biological replicates:

```text
SE = SD / sqrt(n)
```

This follows the convention used in the companion manuscript's figures.

### Stationary-phase detection

The original analysis scripts use two separate parameters to gate independent stationary-phase-like detection for SusI and VI. susiR unifies these into a single `tc`, detected once and reused for both.

This is immaterial when both parameters share their default value of 1.0 h, as used throughout the companion manuscript.

---

## Recent updates

Several capabilities have been added after the initial release:

### Grouped analysis

`host_col` in `run_susi()` and **Analyze in separate groups** in the Shiny application allow a single workbook to contain multiple hosts/strains, each with its own control.

Controls can be auto-detected from common naming patterns or labels containing `"host"`, with a manual override available in the application.

### Visual plate-map editor

The Shiny application provides an interactive alternative to preparing the `name` sheet manually. Users can assign sample, group, and biological-replicate labels directly on a visual 96-well plate layout.

### Raw-format import

`susi_import_raw()` converts exports from different plate readers, including differences in:

- sheet names
- header rows
- well-ID formats
- time units
- non-contiguous well layouts
- English/Japanese headers

### Export customization

The Shiny application allows plot-download dimensions (width/height) and format (PNG/SVG) to be customized, which is useful for multi-panel diagnostic figures and publication preparation.

For implementation details, see:

- `R/run_susi.R`
- `R/raw_import.R`
- `inst/shiny_app/app.R`

---

## Repository structure

```text
susiR/
├── DESCRIPTION
├── R/
│   ├── utils.R              # moving average, integration, well-ID parsing
│   ├── data_input.R         # read_plate_data(), susi_resolve_bio_rep()
│   ├── time_format.R        # Excel time/duration cell detection
│   ├── conditions.R         # condition and MOI handling
│   ├── detection.R          # t0/ti/tc detection and default parameters
│   ├── metrics.R            # SusI, VI, SupI, and ti/tc calculations
│   ├── run_susi.R           # run_susi() -- main entry point
│   ├── diagnose.R           # diagnose_conditions() -- detection QC
│   ├── raw_import.R         # susi_import_raw() -- raw-export conversion
│   ├── plotting.R           # per-condition diagnostic plots
│   ├── plotting_summary.R   # combined and summary plots
│   └── shiny_app.R          # run_app() launcher
├── inst/
│   ├── extdata/
│   │   └── T1_rawdata.xlsx  # bundled example dataset
│   └── shiny_app/
│       ├── app.R
│       └── deploy.R
└── README.md
```

---

## Citation

If you use susiR, please cite:

> Kaneko T, et al. *susiR: standardized, automated quantification of phage lytic sustainability, virulence, and suppression from microplate growth-curve data.* Manuscript in preparation.

---

## License

MIT — see [LICENSE](LICENSE).
