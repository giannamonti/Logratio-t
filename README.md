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
* **Equivalence of logratio transformations** (clr, ilr, alr): all three yield numerically identical Mahalanobis distances and outlier classifications (max difference 
$< 5\times 10^{-14}$).
* **Simulation study** assessing the performance of each method across different sample sizes, dimensions, and contamination levels.

## 🛠 Required Packages

The scripts require the following R libraries:

* **Compositional Data:** `compositions`, `robCompositions`
* **Robust Statistics:** `robustbase`, `fitHeavyTail`, `MASS`, `ContaminatedMixt`
* **Diagnostics & Plotting:** `caret`, `pROC`, `ggplot2`, `patchwork`, `ggVennDiagram`
* **Utilities:** `dplyr`, `tidyr`, `tictoc`, `Surrogate`, `mnormt`, `mvtnorm`, `Rfast`, `ConfidenceEllipse`
* **Dataset:** `StatDA`

You can install all dependencies with:

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
* `exec/utils.R`: Helper functions for confidence ellipses in ternary and coordinate plots.
* `exec/fit_mvt_clr.R`: **Core function.** Implements generalized Student-t fitting for singular clr inputs using Moore-Penrose pseudoinverse.
* `exec/Kola_Pollution.R`: Analysis of the Co-Cu-Ni subset and full 12-part Kola dataset; generates Figures 3 and 4.
* `exec/Kola_outliers.R`: Full-scale outlier detection, consensus analysis, and geochemical interpretation.
* `exec/simulations.R`: Monte Carlo study comparing robust methods across various contamination scenarios.
```

## 🔑 Key Result: Equivalence of Logratio Transformations

A central finding is that CLR, ILR, and ALR are **numerically equivalent** for distribution fitting and outlier detection. Using our generalized `fit_mvt_clr()` function:

* **Mahalanobis Distances:** Max difference $|ilr - clr| < 5 \times 10^{-14}$ (floating-point noise).
* **Outlier Classification:** Symmetric difference between ilr and clr sets is $\emptyset$ (perfect identity).
* **Unified Framework:** The model handles singular (clr) and non-singular (ilr/alr) matrices interchangeably.

## 🚀 Quick Start
1. Install dependencies using the script in the "Required Packages" section.
2. Run `exec/Kola_Pollution.R` to reproduce the main figures and verification tests.
