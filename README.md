# NMA network structure analysis

This standalone R project reconstructs the empirical network sample, simulates a qualified null sample, compares network metrics, compares the sample with the 2000–2019 database, and examines differences by Clinical Topic. Rendering the R Markdown files refreshes their results from the included data.

## Open and run

Open `NMA network structure analysis.Rproj` in RStudio. Install the required R packages if needed:

```r
install.packages(c(
  "rmarkdown", "knitr", "readxl", "igraph", "dplyr", "tidyr",
  "ggplot2", "patchwork", "tnet", "copula", "scales", "stringr", "purrr"
))
```

Render the reports in this order from the project root:

```r
rmarkdown::render("Report 1 - Empirical sample analysis/empirical_database.Rmd")
rmarkdown::render("Report 2 - Null-model simulation and basic characteristics/null_model_sample.Rmd")
rmarkdown::render("Report 3 - Network metric comparisons/metrics comparison.Rmd",
                  output_file = "metrics-comparison.html",
                  output_dir = "Report 3 - Network metric comparisons",
                  output_format = "html_document")
rmarkdown::render("Report 4 - Empirical sample versus old database/empirical_samples_vs_old_database.Rmd")
rmarkdown::render("Report 5 - Clinical topic comparisons/field_class.Rmd")
```

Reports 4 and 5 require the empirical incidence networks produced by Report 1; they can be run independently after Report 1. Each report resolves the project root relative to its own folder, so it can also be rendered from RStudio. Reports 1 and 2 write their reusable RDS and CSV results to the empirical and simulated data folders. Report 3 writes its combined result to its own `results` folder. Report 4 reuses its supplied old-database RDS results by default; set `rebuild_old_database <- TRUE` in its setup chunk to reconstruct them from the included workbooks and extraction table. Report 5 writes its comparison CSVs to its `results` folder. Figures and tables appear in the HTML reports; the R Markdown files do not save separate plot files.

## Project contents

| Folder                                                        | Purpose                                                                                                   | Inputs                                         | Shared functions                               | HTML report               |
| ------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- | ---------------------------------------------- | ---------------------------------------------- | ------------------------- |
| `Empirical data/raw/`                                         | 300 original Excel workbooks used for the empirical sample                                                | —                                              | —                                              | —                         |
| `Empirical data/results/`                                     | Incidence networks and derived empirical results used by Reports 2 and 3                                  | Raw workbooks                                  | —                                              | —                         |
| `Empirical data/old_database/`                                | 303 historical workbooks, extraction table, and cached results for Report 4                                | Historical workbooks and extraction table      | —                                              | —                         |
| `Empirical data/metadata/`                                    | Clinical Topic classification used by Report 5                                                             | Network summary metadata                       | —                                              | —                         |
| `Simulated data/results/`                                     | The published 300-network null sample, the rewiring diagnostic used in Report 2, and derived null results | Empirical results                              | —                                              | —                         |
| `Common R functions/`                                         | Data conversion, network metric computation, and null simulation                                          | —                                              | —                                              | —                         |
| `Report 1 - Empirical sample analysis/`                       | Rebuild trial-by-treatment incidence matrices and describe observed networks                              | Original Excel workbooks                       | `data_processing.R`, `metric_functions.R`      | `empirical_database.html` |
| `Report 2 - Null-model simulation and basic characteristics/` | Fit the joint size model, generate or reuse qualified null networks, and check basic characteristics      | Empirical results; supplied 300-network sample | `simulation_functions.R`, `metric_functions.R` | `null_model_sample.html`  |
| `Report 3 - Network metric comparisons/`                      | Compare empirical and null metrics, fit relationships, and show tables and plots                          | Results from Reports 1 and 2                   | Reads their computed results                   | `metrics-comparison.html` |
| `Report 4 - Empirical sample versus old database/`            | Compare the 2023 empirical sample with the 2000–2019 database                                             | Empirical incidence networks; historical data  | `data_processing.R`, `metric_functions.R`      | `empirical_samples_vs_old_database.html` |
| `Report 5 - Clinical topic comparisons/`                      | Compare empirical network metrics across Clinical Topics                                                  | Empirical incidence networks; topic metadata   | `metric_functions.R`                           | `field_class.html`        |

The incidence matrices have trials as rows, treatments as columns, and binary entries for trial-treatment membership. `data_processing.R` contains the four functions needed from the original external `graph_functions.R`; the source file was not changed. The 300-network RDS cache lets Report 2 reuse the published sample. Removing it causes Report 2 to simulate a new sample with its fixed seed.

The HTML files are provided reference outputs. Render them again after changing code or data. The reports require R and Pandoc in addition to the packages above. The two original R Markdown files in `Null Model/` remain unchanged.
