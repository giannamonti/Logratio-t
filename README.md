# The logratio Student’s t-distribution: a robust model for compositional data analysis

This repository contains the **R code** to reproduce the statistical analyses and results presented in the manuscript:

> **Monti, G. S. et al.** *"The logratio Student's t-distribution: a robust model for compositional data analysis"*.

The project implements the Student's t-distribution within the Logratio approach for Compositional Data (CoDa), providing a robust framework for outlier detection and handling heavy-tailed distributions.

## 📌 Overview
The analysis focuses on:

* **Robust parameter estimation** using the Student's t-distribution.
* **Leave-One-Out (LOO) diagnostics** and Atypicality Index calculation.
* **Comparative analysis** with other robust methods (MCD, COMCoDa, Contaminated Normal).
* **Geochemical interpretation** of consensus outliers (Kola Project dataset).

## 🛠 Required Packages

The script requires the following R libraries:

* **Compositional Data:** `compositions`, `robCompositions`
* **Robust Statistics:** `robustbase`, `fitHeavyTail`, `MASS`, `ContaminatedMixt`
* **Diagnostics & Plotting:** `caret`, `pROC`, `ggplot2`
* **Dataset:** `StatDA`

You can install all dependencies with this command:
```r
install.packages(c("compositions", "robCompositions", "robustbase", "fitHeavyTail", 
                   "MASS", "caret", "pROC", "ggplot2", 
                   "StatDA", "ContaminatedMixt", "dplyr", "tidyr"))
```                  
                   
## 📂 Repository Structure

* `Kola_outliers.R`: Full-scale outlier detection using a consensus of 6 robust methods (Tables 2-3).
* `Kola_Pollution.R`: Detailed analysis of the Co-Cu-Ni subset and PCA/Biplot visualizations (Figures 3-4).
* `utils.R`: Helper functions, including custom confidence ellipse calculations for the Student-t distribution.                   
                   
                   
                   