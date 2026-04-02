set.seed(123)

source("exec/fit_mvt_clr.R") # fit_mvt_clr(), dmt_clr(), maha_clr(),
# pseudo_log_det(), eff_rank()

library(compositions)
library(robustbase)
library(fitHeavyTail)
library(MASS)
library(caret)
library(pROC)
library(ggplot2)
library(patchwork)
library(ggVennDiagram)
library(StatDA)
library(dplyr)

## ---- Load data ----
data(chorizon)
Kola <- chorizon[, c("Co", "Cu", "Ni", "Mg", "Na", "S",
                     "As", "Bi", "Cd", "Sb", "Ag", "Pb")]
dim(Kola)

# =============================================================================
# CLR and ILR are run in parallel to demonstrate numerical equivalence.
# Mahalanobis distances differ by < 5e-14; consensus outlier sets are identical
# (symmetric difference = 0). ALR is not run separately — by transform-
# invariance it yields identical results (see Section 3 of the paper).
# MCD, COMCoDa and CN require full-rank input and are fitted on ILR.
# =============================================================================

Xcoord_ilr <- as.matrix(ilr(Kola))    # N x (D-1), full rank
Xcoord_clr <- as.matrix(clr(Kola))    # N x D,     rank D-1

N    <- nrow(Xcoord_ilr)
d    <- ncol(Xcoord_ilr)   # effective dimension = D-1 = 11
n_df <- N - 1

cat("Dimensions ILR:", dim(Xcoord_ilr), "\n")
cat("Dimensions CLR:", dim(Xcoord_clr), "\n")

# =============================================================================
# 1. Leave-One-Out (LOO) loop — ILR & CLR in parallel
# =============================================================================

d2_loo_norm_ilr <- numeric(N)
d2_loo_t_ilr    <- numeric(N)
atyp_pvals_ilr  <- numeric(N)
nu_loo_ilr      <- numeric(N)

d2_loo_norm_clr <- numeric(N)
d2_loo_t_clr    <- numeric(N)
atyp_pvals_clr  <- numeric(N)
nu_loo_clr      <- numeric(N)

message("Starting LOO loop (ILR + CLR) — ", N, " iterations...")

for (i in 1:N) {
  
  # ---- ILR ------------------------------------------------------------------
  X_train_ilr <- as.matrix(Xcoord_ilr[-i, ])
  z_ilr        <- as.matrix(Xcoord_ilr[i, , drop = FALSE])
  
  fit_n_ilr  <- mvnorm.mle(X_train_ilr)
  mu_n_ilr   <- fit_n_ilr$mu
  sig_n_ilr  <- fit_n_ilr$sigma
  
  qy_ilr <- (1 / (1 + 1/N)) * mahalanobis(z_ilr, mu_n_ilr, sig_n_ilr)
  atyp_pvals_ilr[i]  <- pbeta(qy_ilr / (qy_ilr + n_df), d/2, (n_df - d + 1)/2)
  d2_loo_norm_ilr[i] <- mahalanobis(z_ilr, mu_n_ilr, sig_n_ilr)
  
  fit_t_ilr       <- fit_mvt(X_train_ilr, nu = "iterative", nu_iterative_method = "ECM")
  d2_loo_t_ilr[i] <- mahalanobis(z_ilr, fit_t_ilr$mu, fit_t_ilr$scatter)
  nu_loo_ilr[i]   <- fit_t_ilr$nu
  
  # ---- CLR ------------------------------------------------------------------
  X_train_clr <- Xcoord_clr[-i, ]
  z_clr        <- matrix(Xcoord_clr[i, ], nrow = 1)
  
  mu_n_clr   <- colMeans(X_train_clr)
  sig_n_clr  <- cov(X_train_clr)
  Sinv_n_clr <- ginv(sig_n_clr)
  
  qy_clr <- (1 / (1 + 1/N)) * maha_clr(z_clr, mu_n_clr, Sinv_n_clr)
  atyp_pvals_clr[i]  <- pbeta(qy_clr / (qy_clr + n_df), d/2, (n_df - d + 1)/2)
  d2_loo_norm_clr[i] <- maha_clr(z_clr, mu_n_clr, Sinv_n_clr)
  
  fit_t_clr       <- fit_mvt_clr(X_train_clr)
  Sinv_t_clr      <- ginv(fit_t_clr$scatter)
  d2_loo_t_clr[i] <- maha_clr(z_clr, fit_t_clr$mu, Sinv_t_clr)
  nu_loo_clr[i]   <- fit_t_clr$nu
  
  if (i %% 100 == 0) cat("  Point", i, "/", N, "\n")
}

cat("\nMean nu (LOO) ILR:", round(mean(nu_loo_ilr), 3), "\n")
cat("Mean nu (LOO) CLR:", round(mean(nu_loo_clr), 3), "\n")

cat("\nEquivalence check (Normal LOO distances):\n")
cat("  Max |ILR - CLR|:", max(abs(d2_loo_norm_ilr - d2_loo_norm_clr)), "\n")
cat("Equivalence check (Student-t LOO distances):\n")
cat("  Max |ILR - CLR|:", max(abs(d2_loo_t_ilr - d2_loo_t_clr)), "\n")

# =============================================================================
# 2. Robust Methods (Full Sample)
# =============================================================================

mcd_res    <- covMcd(Xcoord_ilr)
mcd_scores <- mcd_res$mah

comed_res    <- covComed(Xcoord_ilr, reweight = TRUE)
comed_scores <- mahalanobis(Xcoord_ilr, comed_res$center, comed_res$cov)

fit_cn  <- ContaminatedMixt::CNmixt(Xcoord_ilr, G = 1, contamination = TRUE,
                                    model = "VVV", verbose = FALSE)
pred_cn <- factor(ifelse(as.numeric(fit_cn$models[[1]]$v) < 0.5,
                         "Atypical", "Normal"))

# =============================================================================
# 3. Thresholds & Prediction Tables
# =============================================================================
chi2_thresh <- qchisq(0.95, d)

final_t_ilr  <- fit_mvt(Xcoord_ilr, nu = "iterative", nu_iterative_method = "ECM")
gdl_ilr      <- final_t_ilr$nu
thresh_t_ilr <- d * qf(0.95, d, gdl_ilr)

final_t_clr  <- fit_mvt_clr(Xcoord_clr)
gdl_clr      <- final_t_clr$nu
thresh_t_clr <- d * qf(0.95, d, gdl_clr)

cat("\nFull-sample nu  ILR:", round(gdl_ilr, 3),
    "| CLR:", round(gdl_clr, 3), "\n")

preds_ilr <- data.frame(
  ID          = 1:N,
  Atypicality = factor(ifelse(atyp_pvals_ilr >= 0.95,         "Atypical", "Normal")),
  Norm_LOO    = factor(ifelse(d2_loo_norm_ilr >= chi2_thresh,  "Atypical", "Normal")),
  T_LOO       = factor(ifelse(d2_loo_t_ilr    >= thresh_t_ilr, "Atypical", "Normal")),
  MCD         = factor(ifelse(mcd_scores      >  chi2_thresh,  "Atypical", "Normal")),
  COMCoDa     = factor(ifelse(comed_scores    >  chi2_thresh,  "Atypical", "Normal")),
  CN          = pred_cn
)

preds_clr <- data.frame(
  ID          = 1:N,
  Atypicality = factor(ifelse(atyp_pvals_clr >= 0.95,         "Atypical", "Normal")),
  Norm_LOO    = factor(ifelse(d2_loo_norm_clr >= chi2_thresh,  "Atypical", "Normal")),
  T_LOO       = factor(ifelse(d2_loo_t_clr    >= thresh_t_clr, "Atypical", "Normal")),
  MCD         = factor(ifelse(mcd_scores      >  chi2_thresh,  "Atypical", "Normal")),
  COMCoDa     = factor(ifelse(comed_scores    >  chi2_thresh,  "Atypical", "Normal")),
  CN          = pred_cn
)

# =============================================================================
# 4. Outlier Counts & Consensus
# =============================================================================

cat("\n=== ILR — Outlier counts per method ===\n")
print(colSums(preds_ilr[, -1] == "Atypical"))

cat("\n=== CLR — Outlier counts per method ===\n")
print(colSums(preds_clr[, -1] == "Atypical"))

outlier_lists_ilr <- lapply(preds_ilr[, -1], function(col) which(col == "Atypical"))
outlier_lists_clr <- lapply(preds_clr[, -1], function(col) which(col == "Atypical"))

strong_idx_ilr <- Reduce(intersect, outlier_lists_ilr)
strong_idx_clr <- Reduce(intersect, outlier_lists_clr)

cat("\n--- CONSENSUS OUTLIERS ---\n")
cat("ILR: n =", length(strong_idx_ilr),
    "(", round(length(strong_idx_ilr)/N*100, 2), "%)\n")
cat("CLR: n =", length(strong_idx_clr),
    "(", round(length(strong_idx_clr)/N*100, 2), "%)\n")
cat("Symmetric difference (ILR vs CLR):",
    setdiff(union(strong_idx_ilr, strong_idx_clr),
            intersect(strong_idx_ilr, strong_idx_clr)), "\n")

# =============================================================================
# 5. Geochemical Interpretation
# =============================================================================

geom_mean <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))

analyze_outliers <- function(idx, label) {
  Kola_a      <- Kola
  Kola_a$Type <- "Typical"
  Kola_a$Type[idx] <- "Strong Outlier"
  
  comp <- Kola_a %>%
    group_by(Type) %>%
    summarise(across(where(is.numeric),
                     list(Med = median, GMean = geom_mean),
                     .names = "{.col}_{.fn}")) %>%
    tidyr::pivot_longer(-Type,
                        names_to  = c("Element", "Stat"),
                        names_sep = "_") %>%
    tidyr::pivot_wider(names_from = Type, values_from = value) %>%
    mutate(Ratio = `Strong Outlier` / Typical) %>%
    arrange(Stat, desc(Ratio))
  
  cat("\n===", label, "— Most enriched elements (Geometric Mean) ===\n")
  comp %>% filter(Stat == "GMean") %>% head(10) %>% print()
}

analyze_outliers(strong_idx_ilr, "ILR consensus outliers")
analyze_outliers(strong_idx_clr, "CLR consensus outliers")