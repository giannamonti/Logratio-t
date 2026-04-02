# The logratio Student's t-distribution: a robust model for compositional data analysis

This repository contains the **R code** to reproduce the statistical analyses and results presented in the manuscript:

> **Monti, G. S. et al.** *"The logratio Student's t-distribution: a robust model for compositional data analysis"*.

The project implements the Student's t-distribution within the Logratio approach for Compositional Data (CoDa), providing a robust framework for outlier detection and handling heavy-tailed distributions.

## 📌 Overview

The analysis focuses on:

* **Robust parameter estimation** using the Student's t-distribution.
* **Leave-One-Out (LOO) diagnostics** and Atypicality Index calculation.
* **Comparative analysis** with other robust methods (MCD, COMCoDa, Contaminated Normal).
* **Geochemical interpretation** of consensus outliers (Kola Project dataset).
* **Equivalence of logratio transformations** (clr, ilr, alr): all three yield numerically identical Mahalanobis distances and outlier classifications (max difference $< 5×10^{-14}$). The clr is handled directly via the Moore-Penrose pseudoinverse and pseudo-determinant, without requiring any change of basis.
* **Simulation study** assessing the performance of each method across different sample sizes, dimensions, and contamination levels.
* **Envelope goodness-of-fit analysis** confirming the superior fit of the logratio t over the logratio normal.

## 🛠 Required Packages

```r
install.packages(c(
  "compositions", "robCompositions", "robustbase", "fitHeavyTail",
  "MASS", "ContaminatedMixt", "caret", "pROC", "ggplot2",
  "patchwork", "ggVennDiagram", "dplyr", "tidyr", "tictoc",
  "Surrogate", "mnormt", "mvtnorm", "Rfast", "ConfidenceEllipse",
  "StatDA", "ICSOutlier"
))
```

## 📂 Repository Structure

```
exec/
├── utils.R            Helper functions: confidence_ellipse_GM() for
│                      Student-t ellipses in coordinate and ternary plots.
│
├── fit_mvt_clr.R      Generalised Student-t fitting for singular clr input.
│                      Implements fit_mvt_clr(), dmt_clr(), maha_clr(),
│                      pseudo_log_det(), eff_rank() using Moore-Penrose
│                      pseudoinverse and pseudo-log-determinant.
│
├── Kola_Pollution.R   Section 1: Co-Cu-Ni subset (D=3) and Section 2:
│                      full 12-part Kola composition. Fits Normal and
│                      Student-t models for clr, ilr and alr. Generates
│                      Figure 3 (2×2 ternary and coordinate plots) and
│                      Figure 4 (clr-PCA biplot). Includes rotation check
│                      confirming clr SVD projection = ilr in data-driven basis.
│
├── Kola_outliers.R    Full-scale outlier detection on the 12-part Kola dataset.
│                      LOO distances and Atypicality Index for clr and ilr.
│                      Robust methods: MCD, COMCoDa, Contaminated Normal.
│                      Consensus outlier set (42 observations, symmetric
│                      difference clr vs ilr = ∅). Geochemical interpretation.
│
├── envelope_plot.R    Simulation-based envelope analysis (B=999 replications)
│                      for logratio Normal and Student-t goodness of fit.
│                      Generates Figure 5 (2×2 layout, Pollution and Kola).
│
└── simulations.R      Monte Carlo study (100 replications).
                       Scenarios: n ∈ {100,1000}, D ∈ {3,5},
                       contamination ∈ {10%,40%}. Methods compared:
                       T_LOO, N_LOO, Atypicality, MCD, COMCoDa, CN.
                       Metrics: AUC, Sensitivity, Specificity, PPV, NPV.
```


## 🚀 Quick Start

```r
source("exec/fit_mvt_clr.R")   # load generalised fitting functions
source("exec/Kola_Pollution.R") # reproduce Figures 3 and 4
source("exec/Kola_outliers.R")  # reproduce outlier analysis and Table 4
source("exec/envelope_plot.R")  # reproduce Figure 5
```