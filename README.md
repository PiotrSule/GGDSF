# GGDSF: Generalized Gamma Distribution with Shape Function

This repository contains the complete, reproducible R package (`ggdsf`) and source code accompanying the paper:  
**"Generalized Gamma Distribution with Shape Function"** *(Manuscript in preparation / submitted)*.

## Overview

The **Generalized Gamma Distribution with Shape Function (GGDSF)** extends the classical three-parameter Generalized Gamma distribution (Stacy, 1962) by replacing its static power-shape exponent $b$ with a linear time-dependent shape function $w(t) = b + dt$. Under the admissibility condition $b \ge d a e^{-2}$ ($a, b, c > 0$, $d \ge 0$), the GGDSF accommodates L-shaped, unimodal, and bimodal densities as well as flat-bottomed bathtub hazard rate trajectories using only four parameters.

## Installation

You can install the `ggdsf` package directly from GitHub using `pak` (recommended) or `remotes`:

```r
# Recommended installation via pak:
# install.packages("pak")
pak::pak("PiotrSule/GGDSF")

# Alternative via remotes:
# install.packages("remotes")
remotes::install_github("PiotrSule/GGDSF")
```

## Repository Structure

* **`data/aarset.rda`** — Built-in Aarset (1987) electronic device lifetimes dataset ($n = 50$).
* **`R/ggdsf_distribution.R`** — Core distribution functions: CDF (`pggdsf`), survival (`sggdsf`), PDF (`dggdsf`), similarity index (`similarity_ggdsf`), hazard rate (`hggdsf`), reversed hazard rate (`rhggdsf`), quantile (`qggdsf`), random generator (`rggdsf`), and information entropies (`entropy_shannon_ggdsf`, `entropy_renyi_ggdsf`, `entropy_tsallis_ggdsf`).
* **`R/ggdsf_fit.R`** — Maximum likelihood estimation routines for complete (`fit_ggdsf`) and right-censored (`fit_ggdsf_censored`) data with observed Fisher information matrices and LRT diagnostics.
* **`scripts/01_theoretical_shapes.R`** — Reproduces Figures 1–5 (PDF, similarity $M$, HRF, HRAF, and entropy regimes).
* **`scripts/02_monte_carlo_sim.R`** — Reproduces Tables 5–6 (Monte Carlo simulation study evaluating MLE bias, RMSE, empirical/asymptotic SE, and 95% Wald CI coverage).
* **`scripts/03_app_aarset_full.R`** — Section 5.1: Complete Aarset (1987) reliability dataset ($n = 50$, Figure 6).
* **`scripts/04_app_aarset_censored.R`** — Section 5.2: Type-I right-censored Aarset dataset ($C = 75$ weeks, $28.0\%$ censoring, Figure 7).
* **`scripts/05_app_boston_lstat.R`** — Section 5.3: Socioeconomic application to `MASS::Boston$lstat` ($n = 506$, Figure 8).
* **`scripts/06_app_cgd_survival.R`** — Section 5.4: Biomedical survival application to `survival::cgd$tstop` ($n = 203$, Figure 9).

## Requirements

All scripts were developed and tested in **R (>= 4.3.0)** using the following CRAN packages:

```r
install.packages(c("numDeriv", "goftest", "ggplot2", "patchwork", "survival", "MASS"))
```

## Quick Start

```r
library(ggdsf)

# 1. Load the built-in Aarset (1987) dataset (n = 50) and fit GGDSF vs. GGD
data(aarset)
fit_res <- fit_ggdsf(aarset, init = c(80, 1.5, 0.4, 0.06))
print(fit_res$summary_table)
print(fit_res$criteria)
print(fit_res$lrt)

# 2. Evaluate PDF, CDF, Survival, Hazard Rate, Similarity M, and Entropies
t_vals <- c(1, 10, 40, 80)
dggdsf(t_vals, a = 83.3, b = 1.5, c = 0.37, d = 0.06)
pggdsf(t_vals, a = 83.3, b = 1.5, c = 0.37, d = 0.06)
sggdsf(t_vals, a = 83.3, b = 1.5, c = 0.37, d = 0.06)
hggdsf(t_vals, a = 83.3, b = 1.5, c = 0.37, d = 0.06)

similarity_ggdsf(a = 83.3, b = 1.5, c = 0.37, d = 0.06)
entropy_shannon_ggdsf(a = 83.3, b = 1.5, c = 0.37, d = 0.06)

# 3. Generate random variates via numerical quantile inversion
set.seed(42)
sim_sample <- rggdsf(100, a = 83.3, b = 1.5, c = 0.37, d = 0.06)
```

## Citation

If you use this code or the GGDSF distribution in your research, please cite:

> Sulewski, P., Volodin, A. (2026). *Generalized Gamma Distribution with Shape Functions*.
