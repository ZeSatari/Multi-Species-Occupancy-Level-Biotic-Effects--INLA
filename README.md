# Anonymous Code Repository

## Description

This repository contains the R code required to reproduce the analyses presented in the associated manuscript.

## Repository Structure

- `scripts/`: R scripts used in the study.
  - `01_Simulation_detect.R`: Simulates data for the single-species occupancy model with detection-level covariates.
  - `02_Simulation_INLA_detect.R`: Fits the single-species occupancy models using INLA.
  - `03_Simulation_msom.R`: Simulates data for the hierarchical multi-species occupancy model.
  - `04_Simulation_INLA_msom.R`: Fits the hierarchical multi-species occupancy model using INLA.
  - `05_Real_detect.R`: Applies the single-species occupancy model to the empirical dataset.
  - `06_Real_msom.R`: Applies the hierarchical multi-species occupancy model to the empirical dataset.
  - `07_Real_single.R`: Fits separate single-species occupancy models for comparison.

## Empirical Data

The empirical dataset used in the case study is publicly available through the **spOccupancy** R package.

```r
data(hbefTrends, package = "spOccupancy")
```

## Installation

Install the required R packages:

```r
install.packages(c(
  "INLA", "inlabru", "ggplot2", "dplyr", "tidyr",
  "terra", "sf", "tidyverse", "spatstat", "gt",
  "viridis", "viridisLite", "scico", "patchwork",
  "fmesher", "timechange", "lubridate",
  "spOccupancy", "kableExtra"
))
```

Since **INLA** is not available on CRAN, install it from its official repository:

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  )
)
```

## Reproducing the Analysis

Run the scripts in the following order:

```r
# Simulation study: detection-level effects
source("scripts/01_Simulation_detect.R")
source("scripts/02_Simulation_INLA_detect.R")

# Simulation study: hierarchical MSOM
source("scripts/03_Simulation_msom.R")
source("scripts/04_Simulation_INLA_msom.R")

# Empirical analysis: single-species model
source("scripts/05_Real_detect.R")

# Empirical analysis: multi-species model
source("scripts/06_Real_msom.R")

# Separate single-species analyses
source("scripts/07_Real_single.R")
```

## License

This project is released under the MIT License.

Revise README
