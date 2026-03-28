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
* **Equivalence of logratio transformations** (CLR, ILR, ALR): all three yield numerically identical Mahalanobis distances and outlier classifications (max difference < 5×10$^{-14}$).
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
exec/
├── utils.R              # Helper functions: confidence_ellipse_GM() for
│                        # Student-t ellipses in coordinate and ternary plots.
│
├── fit_mvt_clr.R        # Generalised Student-t fitting for singular CLR input.
│                        # Implements fit_mvt_clr(), dmt_clr(), maha_clr(),
│                        # pseudo_log_det() and eff_rank() using Moore-Penrose
│                        # pseudoinverse (MASS::ginv) and pseudo-log-determinant.
│                        # Handles all logratio transformations (CLR, ILR, ALR)
│                        # in a unified framework without requiring an explicit
│                        # change of basis.
│
├── Kola_Pollution.R     # Section 1: Co-Cu-Ni subset (D=3).
│                        #   - ILR, ALR and CLR fits (Normal and Student-t).
│                        #   - AIC comparison across transformations.
│                        #   - Numerical verification of transformation equivalence.
│                        #   - Figure 3: 2×2 layout (ternary and coordinate plots
│                        #     for ILR and CLR).
│                        # Section 2: Full 12-part Kola composition.
│                        #   - AIC summary table (CLR, ILR, ALR).
│                        #   - Figure 4: CLR-PCA biplot with Normal and Student-t
│                        #     confidence ellipses.
│
├── Kola_outliers.R      # Full-scale outlier detection on the 12-part Kola dataset.
│                        #   - LOO distances and Atypicality Index (ILR and CLR).
│                        #   - Robust methods: MCD, COMCoDa, Contaminated Normal.
│                        #   - Consensus outlier set (all 6 methods agree).
│                        #   - Geochemical interpretation of consensus outliers.
│                        #   - Equivalence check: ILR and CLR yield identical
│                        #     consensus outliers (symmetric difference = ∅).
│
└── Simulations.R        # Monte Carlo simulation study (100 replications).
                         #   - Scenarios: n ∈ {100, 1000}, p ∈ {3, 5},
                         #     contamination ∈ {10%, 40%}.
                         #   - Methods compared: Student-t LOO, Normal LOO,
                         #     Atypicality Index, MCD, COMCoDa, CN.
                         #   - Metrics: AUC, Sensitivity, Specificity, PPV, NPV.
                         #   - ILR coordinates used throughout (CLR and ALR yield
                         #     identical results — see Kola_outliers.R).
```

## 🔑 Key Result: Equivalence of Logratio Transformations

A central finding of this work is that CLR, ILR and ALR are **numerically equivalent** for all analyses performed:

```r
# Normal Mahalanobis distances
Max |ILR - CLR|: 4.97e-14   # floating-point noise only
Max |ILR - ALR|: 4.62e-14

# Consensus outlier sets
Symmetric difference (ILR vs CLR): ∅   # identical sets
```

The `fit_mvt_clr()` function in `exec/fit_mvt_clr.R` implements a generalised version of `fitHeavyTail::fit_mvt()` that works directly on the full D-dimensional singular CLR matrix, using the Moore-Penrose pseudoinverse and pseudo-log-determinant. This avoids any explicit change of basis and provides a unified treatment for all logratio transformations.