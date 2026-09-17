# susiR

**susiR** computes four complementary indices of bacteriophage lytic activity from microplate (optical density) growth-curve time-course data:

- **SusI** (Sustainability Index) — duration and depth of lysis suppression, from lysis onset (t0) to resistance emergence (ti)
- **VI** (Virulence Index) — early infection dynamics, following Storms et al. (2020)
- **SupI** (Suppression Index) — overall growth suppression over a defined observation window, following Kim et al. (2024)
- **ti/tc** — a time-ratio metric relating suppression duration to the control culture's normal growth cycle

susiR is the companion software for the Sustainability Index manuscript (Kaneko et al., in preparation). Unlike the original analysis scripts developed for that study, susiR makes no assumptions about experimental design: the number of conditions, their labels, and the number and layout of biological replicates are all inferred from a user-supplied sample-mapping table at run time, so the same code handles arbitrary plate layouts without modification.

## Key features

- **Flexible replicate count.** Biological replicates are not fixed at 3. An explicit `bio_rep` column in the mapping sheet is used exactly as given — 1, 2, 4, or any other number, arranged however the wells were laid out.
- **Flexible replicate layout.** Because replicate membership can be declared per well, replicates need not occupy a fixed, contiguous block of plate columns. Without a `bio_rep` column, susiR falls back to a positional convention (consecutive blocks of 12 wells per condition, in plate order), adapting automatically to however many wells a condition has.
- **Flexible condition labels.** Conditions are never assumed to be an MOI dilution series. A single sheet can equally contain several phages tested at one fixed MOI, or any other labeling scheme; every distinct label in the mapping sheet is treated as one condition. The MOI-dependent global VI/MV50 dose-response fit is computed only when condition labels parse as a numeric MOI series (checked automatically); local, per-condition SusI/VI/SupI/ti-tc are always computed.
- **Import from any plate reader.** Raw exports in unfamiliar layouts (unknown sheet name, header row, well-ID spelling, time units) can be auto-converted into the format susiR expects — see "Importing from other plate readers" below.
- **Point-and-click interface.** A bundled Shiny application provides the full workflow — upload, automatic condition/replicate detection, parameter tuning, diagnostic plots, and results download — without writing R code.

## Getting started

No package installation step is required to try susiR — source the files in `R/` directly:

```r
for (f in list.files("R", full.names = TRUE)) source(f)
```

(A standard package build via `devtools::load_all()` or `R CMD build` also works, since the functions are documented with standard roxygen-style comments.)

## Input format

**OD sheet** (default name `"R"`): first column `Time`, every other column one well. `Time` may be a plain number (interpreted as elapsed seconds) or a genuine Excel time/duration cell (e.g. formatted as `[h]:mm`); `read_plate_data()` inspects the workbook's actual cell formatting to distinguish these automatically (`time_unit = "auto"`, the default), since `openxlsx` returns a different raw number for each and treating one as the other would silently corrupt every downstream time value by a factor of 86400.

**Mapping sheet** (default name `"name"`): one row per well.

| Column | Required? | Meaning |
|---|---|---|
| `num` | yes | Well ID (must match a column in the OD sheet) |
| `sample` | yes | Condition label — any text, no fixed vocabulary |
| `bio_rep` | no | Explicit biological-replicate ID. Omit entirely to use the positional fallback (blocks of 12 wells, plate order). If present, must be filled in for every well of a given condition. |

Column names are configurable via `well_col`, `condition_col`, `bio_rep_col` if a sheet uses different headers.

### Importing from other plate readers

Raw exports from plate-reader software rarely match the shape above directly — each has its own sheet name, a header row buried under instrument metadata, well IDs spelled differently (`A1`, `A:1`, `A01`, or embedded in a longer label such as `"Sample0001 (B02)"`), and a time column with its own units. `susi_import_raw()` converts an arbitrary raw export into the expected shape:

```r
susi_import_raw("my_raw_export.xlsx", output_path = "converted.xlsx")
#> Imported from sheet 'Results Table' (header row 7): 96 wells, time
#> column 'Row time (sec)' (unit: seconds, detected via header_text),
#> 0.00 to 24.50 hours over 99 points.
#> Wrote converted.xlsx -- fill in the 'sample' (and optionally 'bio_rep')
#> column in its 'name' sheet, then use it with read_plate_data()/run_susi()
#> as normal.
```

It scans every sheet for a row containing a time label (English or Japanese: "time"/"時間", "sec"/"秒", "min"/"分") and at least six well-ID-like column headers, drops everything else (temperature columns, read-type/read-number columns, trailing metadata), and writes a workbook with a ready-to-use `"R"` sheet plus a `"name"` sheet template (well IDs filled in; `sample`/`bio_rep` left blank, since experimental design assignment always requires a human). The printed summary should be checked against expectations before the result is used further; a low-confidence match fails with an explicit error rather than guessing.

Validated against exports from three different instruments: a BioTek/Agilent Epoch2 reader (multi-row metadata preamble, Excel-duration time column), a Promega reader (distinct column layout, `A:1`-style well IDs, embedded control characters in headers — parsed via the `readxl` package as a fallback where `openxlsx` fails), and a Thermo Fisher SkanIt export in Japanese with a sparse, non-contiguous well selection and well IDs embedded inside longer sample labels. If a file's layout is too unusual for automatic detection, either pass `sheet =` to indicate the correct tab, or convert the file by hand to match the bundled example (`inst/extdata/T1_rawdata.xlsx`).

The Shiny application (below) includes an equivalent "Import a raw plate-reader export" section.

## Usage

```r
for (f in list.files("R", full.names = TRUE)) source(f)

## Bundled example data, no bio_rep column -> positional fallback (3 x 12)
res <- run_susi("inst/extdata/T1_rawdata.xlsx",
                 control_label = "Ct", exclude = "free",
                 calculation_method = "biological_replicates",
                 time_limit_hours = 24)
res$summary        # SusI / VI_local / SupI / ti_tc per condition (+ per bio_rep rows)
res$global_vi      # VI / MV50 across the MOI series (auto-detected as applicable here)

## Custom data with an explicit, irregular replicate structure
## (e.g. 2 biological replicates of different sizes) -- add a `bio_rep`
## column to the mapping sheet; no code changes needed.
res2 <- run_susi("my_plate.xlsx", control_label = "Ct")

## Several phages at one fixed MOI, in a single sheet -- conditions are
## treated as labels; global VI/MV50 (which needs a numeric MOI series) is
## skipped automatically, with an explanatory message, while
## SusI/VI/SupI/ti-tc are still computed per phage.
res3 <- run_susi("cocktail_screen.xlsx", control_label = "Ct")
res3$global_vi$note

## Fast parameter tuning without computing indices:
diag <- diagnose_conditions("inst/extdata/T1_rawdata.xlsx", control_label = "Ct", exclude = "free",
                             output_excel = "T1_diagnosis.xlsx")
```

## Plotting

Functions in `R/plotting.R` and `R/plotting_summary.R` (require `ggplot2` and `patchwork`) turn results into figures:

```r
## Per-condition diagnostic view: control vs. treated OD curves with
## t0/ti/tc marked and the SusI numerator area shaded
plot_condition("inst/extdata/T1_rawdata.xlsx", condition = "10^0",
                control_label = "Ct", time_limit_hours = 24)
plot_all_conditions("inst/extdata/T1_rawdata.xlsx", control_label = "Ct", exclude = "free",
                     time_limit_hours = 24, output_file = "diagnostic_overview.png")

## Every condition's mean curve on one plot, as mean +/- SD across
## biological replicates
plot_combined_curves("inst/extdata/T1_rawdata.xlsx", control_label = "Ct", exclude = "free",
                      time_limit_hours = 24)

## The indices as charts (2x2 grid: SusI/VI/SupI/ti-tc)
res <- run_susi("inst/extdata/T1_rawdata.xlsx", control_label = "Ct", exclude = "free",
                 calculation_method = "biological_replicates", time_limit_hours = 24)
plot_metrics_summary(res)

## Well- and replicate-level view: individual wells (circles), each
## biological replicate's mean (coloured triangle), and the condition mean
## +/- SE (black diamond and error bar) -- within- and between-replicate
## variability visible together, in the same style used in the companion
## manuscript's figures.
plot_metric_superplot("inst/extdata/T1_rawdata.xlsx", control_label = "Ct", exclude = "free",
                       time_limit_hours = 24)
plot_metric_superplot("inst/extdata/T1_rawdata.xlsx", metrics = "SusI",
                       control_label = "Ct", exclude = "free", time_limit_hours = 24)
```

`plot_combined_curves()`, `plot_metrics_summary()`, and `plot_metric_superplot()` adapt automatically to the condition labels: a numeric MOI series is colour-graded and positioned by log10(MOI); any other labeling scheme falls back to categorical colour/bar/point charts. `plot_metric_superplot()` runs `run_susi()` internally at both the `individual_wells` and `biological_replicates` levels, so it takes somewhat longer than the other plotting functions.

## Shiny application

A point-and-click interface is provided at `inst/shiny_app/app.R`, covering the full workflow — upload, automatic condition/replicate detection, interactive parameter tuning with live diagnostic plots, results table, and Excel/PNG downloads.

Requires three additional packages: `install.packages(c("shiny", "DT", "patchwork"))`.

**To launch:**
- Open `inst/shiny_app/app.R` in RStudio and click the **Run App** button that appears in the editor pane; or
- From the package root: `source("R/shiny_app.R"); run_app()`

**Workflow:** upload a `.xlsx` file → the control/exclude dropdowns populate automatically from the mapping sheet → adjust detection parameters in the sidebar if needed → **Run analysis** → review the *Overview*, *Combined view*, and *Metrics summary* tabs (plus *Condition detail* for anything that looks off) → download the results table and any figures.

### Deployment

Running the app locally requires R installed. For unrestricted browser access with no local install, deploy to [Posit Connect Cloud](https://connect.posit.cloud) — see `inst/shiny_app/deploy.R` for setup instructions (either via the `rsconnect` R package, or by publishing directly from a GitHub repository).

## Validation

SusI and VI have been checked against the exact values reported in the companion manuscript for the T1/*E. coli* MG1655 dataset at MOI = 1:

| Index | susiR | Manuscript |
|---|---|---|
| SusI | 0.664 | 0.664 |
| VI (local, MOI = 1) | 0.86 | 0.84 |
| VI (global) | 0.84 | 0.84 |

ti/tc and SupI were verified by re-deriving their calculation formulas directly from the original analysis scripts rather than against a single reported value; both are implemented to match exactly. The remaining five MOI conditions in the dataset have not been checked against exact reported numbers — the manuscript reports precise values only at MOI = 1, with the full dilution series shown graphically (Fig. 2) rather than tabulated — but a comparison against that figure shows close agreement across all four indices and all MOI conditions for the T1 series.

### Implementation notes

A few calculation details are worth stating explicitly for reproducibility:

- **Control stationary-phase time (tc)** is computed once per dataset, from all control wells pooled together, and reused throughout: as the denominator bound for SusI, the integration bound for VI, and the fallback endpoint for ti. SusI's numerator (the t0→ti window) uses each biological replicate's own control curve; only the denominator uses the pooled curve.
- **VI** is integrated from t = 0 to the same control-stationary-phase time (tc) used for SusI, following Storms et al. (2020), rather than over a fixed window to the end of the observation period.
- **SupI's** integration window defaults to `time_limit_hours` rather than a fixed 30 hours, matching how the original pipeline was run for the companion manuscript.
- **Reported variability** (`SusI_se`, `VI_local_se`, `SupI_se`, `ti_tc_se` in `run_susi()$summary`, and the error bars in `plot_metric_superplot()`) is the standard error across biological replicates (SD / √n), matching the convention used in the companion manuscript's figures.
- The original analysis scripts use two separate (usually identically-valued) parameters to gate independent stationary-phase-like detection for SusI and VI; susiR unifies these into a single `tc`, detected once and reused for both. This is immaterial when both parameters share their default value (1.0 h, as used throughout the companion manuscript).

## Repository structure

```
susiR/
├── DESCRIPTION
├── R/
│   ├── utils.R              # moving average, trapezoidal integration, well-ID parsing
│   ├── data_input.R         # read_plate_data(), susi_resolve_bio_rep()
│   ├── time_format.R        # Excel time/duration cell detection
│   ├── conditions.R         # susi_resolve_conditions(), susi_check_moi_applicable()
│   ├── detection.R          # susi_detect_t0/ti/tc(), susi_default_params()
│   ├── metrics.R            # susi_calc_susi/vi_local/global_vi/supi/time_ratio()
│   ├── run_susi.R           # run_susi() -- main entry point
│   ├── diagnose.R           # diagnose_conditions() -- per-well detection QC
│   ├── raw_import.R         # susi_import_raw() -- convert arbitrary plate-reader exports
│   ├── plotting.R           # per-condition diagnostic plots
│   ├── plotting_summary.R   # combined-curve and per-metric summary plots
│   └── shiny_app.R          # run_app() launcher
├── inst/
│   ├── extdata/T1_rawdata.xlsx   # bundled example dataset
│   └── shiny_app/                # Shiny application (app.R, deploy.R)
└── README.md
```

## Citation

If you use susiR, please cite:

> Kaneko T, et al. susiR: standardized, automated quantification of phage lytic sustainability, virulence, and suppression from microplate growth-curve data. *Manuscript in preparation.*

## License

MIT — see [LICENSE](LICENSE).
