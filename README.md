# Multi-Species Occupancy Models with Detection-Level Biotic Effects - INLA

## Description
This repository contains the complete code to reproduce all analyses, simulations, and figures from the paper:  
**"Spatio-Temporal Multi-Species Occupancy Models with Detection-Level Biotic Effects Using INLA"**.

## Repository Structure

-   `scripts/`: All R scripts, organized by analysis:
    -   `01_Simulation_detect.R`: Generates synthetic data for the single-species occupancy model.
    -   `02_Simulation_INLA_detect.R`: Runs the single-species occupancy model comparisons with detection-level biotic effects.
    -   `03_Simulation_msom.R`: Generates synthetic data for the hierarchical Multi-Species Occupancy Model (MSOM).
    -   `04_Simulation_INLA_msom.R`: Runs the hierarchical Multi-Species Occupancy Model comparisons.
    -   `05_Real_detect.R`: Fits the single-species occupancy model with detection-level biotic effects.
    -   `06_Real_msom.R`: Fits the hierarchical Multi-Species Occupancy Model (MSOM) to real data.
    -   `07_Real_single.R`: Fits the single-species occupancy models to real data.

### Empirical Data (Case Study)
The real-world data used in the case study (Hubbard Brook bird surveys) were sourced from the `spOccupancy` R package.
 # Load the dataset
    data(hbefTrends, package = "spOccupancy")

## Installation & Setup

### 1. Install Required R Packages
Run the following command in your R session to install all necessary dependencies:

```r
install.packages(c("INLA", "inlabru", "ggplot2", "dplyr", "tidyr", "terra", "sf", "tidyverse",
                   "spatstat", "gt", "viridis", "viridisLite", "scico", "patchwork", "fmesher",
                   "timechange", "lubridate", "spOccupancy", "kableExtra"))

Important Note: The INLA package is not on CRAN. Install it from its dedicated repository:
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))


2. Reproduce the Analysis
To replicate the main findings and figures from the paper, execute the scripts in order:

r
# 1. Run the Single-Species simulation study: Detection-based biotic information
source("scripts/01_Simulation_detect.R")
source("scripts/02_Simulation_INLA_detect.R")


# 2. Run the Multi-Species Simulation Study
source("scripts/03_Simulation_msom.R")
source("scripts/04_Simulation_INLA_msom.R")


# 3. Run the real-world case study: Detection-based biotic information
source("scripts/05_Real_detect.R")


# 4. Run the real-world case study: Multi-Species Occupancy Model
source("scripts/06_Real_msom.R")
source("scripts/07_Real_single.R")

##Citation
This code supports the research presented in our manuscript on spatio-temporal multi-species occupancy modeling. The manuscript is currently under peer review.

For now, if you use or adapt this code, please:

1- Acknowledge this GitHub repository

2- Provide a link to this repo in your work

3- Cite the associated publication once it becomes available

##Contact
For questions or issues regarding the code, please open an Issue on this GitHub repository or contact the corresponding author at: ze.satari@gmail.com.

##License
This project's code is released under the MIT License.




