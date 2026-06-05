# Artemis-II-Data-Methodology
Bayesian within-subject longitudinal modeling applied to NHANES proxy data to demonstrate a methodology for extracting physiological insights from small-sample, high-dimensional spaceflight data. Submitted to NASA for Artemis II Human Research Data Methodology.

## Component 2: Methodology Demonstration — README

---

## Overview

This repository contains the full demonstration of a **Bayesian within-subject longitudinal modeling methodology** applied to NHANES proxy data.

The core scientific problem mirrored here is extracting meaningful physiological trajectory estimates from a very small, high-performing cohort (n=4 astronauts) with multiple biomarker measurements across pre-flight, in-flight, and post-flight timepoints. The methodology uses informative priors derived from a large terrestrial baseline population (NHANES) to compensate for the inherent statistical limitations of a four-person dataset.

**All data is downloaded automatically from publicly available CDC servers. No registration, account, or manual data download is required.**

---

## Repository Structure

```
artemis-ii-methodology/
│
├── analysis.R              # Main analysis script (run this)
├── README.md               # This file
│
└── outputs/                # Auto-generated when script runs
    ├── 01_cohort_summary.csv
    ├── 02_prior_distributions.png
    ├── 03_trajectory_plot.png
    ├── 04_posterior_summary.csv
    ├── 05_credible_intervals.png
    └── 06_sensitivity_analysis.png
```

---

## System Requirements

| Component | Requirement |
|-----------|-------------|
| R version | 4.2.0 or higher (4.3.x recommended) |
| Operating System | Windows 10/11, macOS 12+, or Linux (Ubuntu 20.04+) |
| RAM | Minimum 8 GB (16 GB recommended for MCMC sampling) |
| Disk space | ~500 MB (for packages and outputs) |
| Internet connection | Required (for NHANES data download and package installation) |
| Estimated runtime | 5–15 minutes depending on hardware |

---

## Dependencies

All R packages are installed automatically by the script if not already present. The following packages are required:

| Package | Purpose |
|---------|---------|
| `nhanesA` | NHANES data download from CDC |
| `brms` | Bayesian regression models via Stan |
| `tidyverse` | Data manipulation and visualization |
| `bayesplot` | Bayesian model visualization |
| `ggplot2` | Plotting |
| `dplyr` | Data wrangling |
| `tidyr` | Data reshaping |
| `patchwork` | Plot composition |
| `posterior` | Posterior summaries |
| `loo` | Leave-one-out cross-validation |

> **Note:** `brms` depends on **Stan** (via the `rstan` or `cmdstanr` backend). If you have not used `brms` before, the first run may prompt you to install the Stan toolchain. See the Troubleshooting section below for guidance.

---

## Step-by-Step Reproduction Instructions

### Step 1: Install R

Download and install R (version 4.2.0 or higher) from the official CRAN mirror:  
https://cran.r-project.org/

RStudio is recommended but not required:  
https://posit.co/download/rstudio-desktop/

---

### Step 2: Clone or Download This Repository

**Option A — GitHub (recommended):**
```bash
git clone https://github.com/[your-username]/artemis-ii-methodology.git
cd artemis-ii-methodology
```

**Option B — Manual download:**  
Download the ZIP from the repository page, extract it, and open the folder.

---

### Step 3: Open the Script

Open `analysis.R` in RStudio or any R-compatible editor.

Confirm your working directory is set to the repository root (the folder containing `analysis.R`). In RStudio, go to:  
**Session → Set Working Directory → To Source File Location**

Or in the R console:
```r
setwd("/path/to/artemis-ii-methodology")
```

---

### Step 4: Run the Script

**Option A — RStudio:**  
Click **Source** (top right of the editor pane) to run the entire script end-to-end.

**Option B — R console:**
```r
source("analysis.R")
```

**Option C — Command line / terminal:**
```bash
Rscript analysis.R
```

The script will:
1. Automatically install any missing packages
2. Download NHANES 2017–2018 data directly from CDC public servers
3. Filter and construct the astronaut-like cohort
4. Derive informative Bayesian priors from the baseline population
5. Construct the small-cohort pseudo-longitudinal dataset
6. Fit three Bayesian mixed-effects models via MCMC (this step takes the longest)
7. Generate and save all output files to the `outputs/` folder

Progress messages are printed to the console at each step.

---

### Step 5: Verify Outputs

After the script completes, confirm the following six files exist in the `outputs/` directory:

| File | Description |
|------|-------------|
| `01_cohort_summary.csv` | Descriptive statistics for the astronaut-like NHANES cohort |
| `02_prior_distributions.png` | Histograms of prior distributions derived from NHANES baseline |
| `03_trajectory_plot.png` | Individual biomarker trajectories across pre-, in-, and post-mission timepoints |
| `04_posterior_summary.csv` | Posterior fixed-effect estimates with 95% credible intervals and probability of direction |
| `05_credible_intervals.png` | Posterior distribution plots for WBC, lymphocyte %, and log(hs-CRP) |
| `06_sensitivity_analysis.png` | Comparison of informative vs. weakly informative prior specifications |

---

## Reproducibility Notes

- A fixed random seed (`set.seed(20260405)`) is set at the beginning of the script and passed to all `brms` model calls (`seed = 20260405`). Results will be numerically identical across runs on the same machine and R version.
- Minor numerical differences in posterior estimates (~0.001) may occur across different operating systems or R versions due to differences in floating-point precision and Stan compiler behavior. These differences do not affect scientific conclusions.
- NHANES data is downloaded directly from CDC public servers at runtime. CDC data files are versioned and stable; the 2017–2018 cycle used here (`DEMO_J`, `CBC_J`, `HSCRP_J`, `PAQ_J`, `SMQ_J`) has been publicly available since 2020 and is not subject to change.
- The `outputs/` directory is created automatically if it does not exist. Do not manually create it beforehand.

---

## Proxy Dataset Information

| Dataset | Source | Variables Used | Access |
|---------|--------|---------------|--------|
| NHANES 2017–2018 Demographics (DEMO_J) | CDC | Age, sex, education, income | Public, no registration |
| NHANES 2017–2018 Complete Blood Count (CBC_J) | CDC | WBC, lymphocyte %, neutrophil %, monocyte %, hemoglobin | Public, no registration |
| NHANES 2017–2018 hs-CRP (HSCRP_J) | CDC | High-sensitivity C-reactive protein | Public, no registration |
| NHANES 2017–2018 Physical Activity (PAQ_J) | CDC | Vigorous/moderate activity indicators | Public, no registration |
| NHANES 2017–2018 Smoking (SMQ_J) | CDC | Lifetime smoking status | Public, no registration |

All data is accessed via the `nhanesA` R package, which interfaces directly with the CDC's public NHANES data portal.

---

## Methodological Summary

The analysis proceeds in six stages:

1. **Data acquisition:** NHANES 2017–2018 physiological and behavioral data downloaded from CDC.
2. **Cohort filtering:** Dataset restricted to an astronaut-like cohort (age 25–55, non-smoker, physically active, college-educated, no acute inflammation) to serve as the terrestrial baseline.
3. **Prior derivation:** Descriptive statistics from the full astronaut-like cohort are used to construct informative Bayesian prior distributions, encoding healthy human physiology into the model before observing the small mission cohort.
4. **Small-cohort construction:** A pseudo-longitudinal dataset (n=6, three timepoints) is constructed from the baseline population, with spaceflight-consistent perturbation magnitudes drawn from published ISS immune research (Crucian et al. 2015; Strangman et al. 2014).
5. **Bayesian model fitting:** Three Bayesian mixed-effects models (one per biomarker: WBC, lymphocyte %, log hs-CRP) are fit using `brms`, each with fixed timepoint effects and by-subject random intercepts and slopes.
6. **Posterior inference:** Posterior distributions, 95% credible intervals, and probability-of-direction statistics are extracted and visualized. A sensitivity analysis compares informative vs. weakly informative priors to confirm prior specifications do not inappropriately dominate the likelihood.

---

## Troubleshooting

**Stan / C++ toolchain not found:**  
`brms` requires a working C++ compiler to interface with Stan. If you receive a toolchain error on first use:
- **Windows:** Install Rtools from https://cran.r-project.org/bin/windows/Rtools/
- **macOS:** Install Xcode Command Line Tools: `xcode-select --install`
- **Linux:** Install `build-essential`: `sudo apt-get install build-essential`

Then restart R and rerun the script.

**NHANES download fails:**  
If the CDC server is temporarily unavailable, wait a few minutes and rerun. The `nhanesA` package downloads directly from CDC public servers and requires an active internet connection.

**MCMC takes longer than expected:**  
The `brms` models use 4 chains with 2000 iterations each. Runtime depends on available CPU cores. The script automatically detects and uses all available cores (`cores = min(4, parallel::detectCores())`). On a 2-core machine, expect runtime closer to 15–20 minutes.

**Package installation errors:**  
If automatic installation fails, install packages manually:
```r
install.packages(c("nhanesA", "brms", "tidyverse", "bayesplot",
                   "patchwork", "posterior", "loo"),
                 repos = "https://cran.rstudio.com/",
                 dependencies = TRUE)
```

---

## Contact

For questions about this submission, contact: Soumya at s.soumyaxd@gmail.com
