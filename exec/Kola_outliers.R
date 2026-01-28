library(compositions)
library(robustbase)
library(fitHeavyTail)
library(MASS)
library(caret)
library(pROC)
library(ggplot2)
library(patchwork)
library(ggVennDiagram)

## ---- Load data ----
data(chorizon) # StatDA package
Kola <- chorizon[, c("Co", "Cu", "Ni", "Mg", "Na", "S",
                     "As", "Bi", "Cd", "Sb", "Ag", "Pb")]
dim(Kola)

# --- 1. Initialization ---
Xcoord <- ilr(Kola)
N <- nrow(Xcoord)
d <- ncol(Xcoord)
n_df <- N - 1

# Vectors to store Leave-One-Out (LOO) distances
d2_loo_norm <- numeric(N)
d2_loo_t    <- numeric(N)
atyp_pvals  <- numeric(N)


for (i in 1:N) {
  X_menoi <- as.matrix(Xcoord[-i, ])
  z <- Xcoord[i, , drop = FALSE]
  
  # A. Normal LOO (for Atypicality Index and Normal LOO Distance)
  fit_norm_i <- mvnorm.mle(X_menoi)
  mu_n <- fit_norm_i$mu
  sig_n <- fit_norm_i$sigma
  
  # Distance for Atypicality (qy formula)
  qy <- (1/(1 + 1/N)) * mahalanobis(z, mu_n, sig_n)
  atyp_pvals[i] <- pbeta(qy / (qy + n_df), d/2, (n_df - d + 1)/2)
  
  # Distance for Normal LOO Confusion Matrix
  d2_loo_norm[i] <- mahalanobis(z, mu_n, sig_n)
  
  # B. Student-t LOO
  fit_t_i <- fit_mvt(X_menoi, nu = "iterative", nu_iterative_method = "ECM")
  d2_loo_t[i] <- mahalanobis(z, fit_t_i$mu, fit_t_i$scatter)
  
  if(i %% 100 == 0) cat("Point", i, "of", N, "\n")
}

# --- 2. Robust Methods (Full Sample - LOO not required here) ---
mcd_res    <- covMcd(Xcoord)
mcd_scores <- mcd_res$mah

comed_res    <- covComed(Xcoord, reweight = TRUE)
comed_scores <- mahalanobis(Xcoord, comed_res$center, comed_res$cov)

fit_cn <- ContaminatedMixt::CNmixt(Xcoord, G = 1, contamination = TRUE, model = "VVV", verbose = FALSE)
pred_cn <- factor(ifelse(as.numeric(fit_cn$models[[1]]$v) < 0.5, "Atypical", "Normal"))


# --- 3. Predictions and Thresholds ---

# T threshold (using the Full model for consistency)
final_t_fit <- fit_mvt(Xcoord, nu = "iterative", nu_iterative_method = "ECM")
gdl_t <- final_t_fit$nu
thresh_t <- d * qf(0.95, d, gdl_t)
t_scores <- mahalanobis(Xcoord, final_t_fit$mu, final_t_fit$scatter)


nu_loo_values <- numeric(N)
# Fast loop to extract nu (without recalculating everything)
message("Extracting nu values via LOO...")
for (i in 1:N) {
  X_menoi <- as.matrix(Xcoord[-i, ])
  # Fast extraction of the nu parameter
  fit_i <- fit_mvt(X_menoi, nu = "iterative", nu_iterative_method = "ECM")
  nu_loo_values[i] <- fit_i$nu
}

mean(nu_loo_values)

final_norm_fit <- mvnorm.mle(as.matrix(Xcoord))
norm_scores <- mahalanobis(Xcoord, center = final_norm_fit$mu, 
                           cov = final_norm_fit$sigma)


preds <- data.frame(
  ID = 1:N,
  Atypicality = factor(ifelse(atyp_pvals >= 0.95, "Atypical", "Normal")),
  Norm_LOO    = factor(ifelse(d2_loo_norm >= qchisq(0.95, d), "Atypical", "Normal")),
  T_LOO       = factor(ifelse(d2_loo_t >= thresh_t, "Atypical", "Normal")),
  MCD         = factor(ifelse(mcd_scores > qchisq(0.95, d), "Atypical", "Normal")),
  COMCoDa     = factor(ifelse(comed_scores > qchisq(0.95, d), "Atypical", "Normal")),
  CN          = pred_cn
)


# --- 4. Validation and Consensus ---
print("Outlier Count (With recalculated T-LOO):")
print(colSums(preds[,-1] == "Atypical"))

# Strong Outliers (Unanimity among T_LOO, MCD, CN, COMCoDa, Norm, and Atypicality)
strong_idx <- which(preds$T_LOO == "Atypical" & 
                      preds$MCD == "Atypical" & 
                      preds$CN == "Atypical" &
                      preds$COMCoDa == "Atypical" &
                      preds$Atypicality == "Atypical" & 
                      preds$Norm_LOO == "Atypical")

cat("Outliers identified by consensus (including T-LOO):", length(strong_idx), "\n")

# --- 2. Results Summary ---

results_df <- data.frame(
  ID = 1:N,
  MCD = preds$MCD,
  COMCoDa = preds$COMCoDa,
  T_Student = preds$T_LOO,
  CN = preds$CN,
  Atypicality = preds$Atypicality,
  Norm =  preds$Norm_LOO
)

print("Count of Identified Outliers:")
print(colSums(results_df[, -1] == "Atypical", na.rm = TRUE))



outlier_lists <- list(
  MCD = which(results_df$MCD == "Atypical"),
  COMCoDa = which(results_df$COMCoDa == "Atypical"),
  T_Student = which(results_df$T_Student == "Atypical"),
  CN = which(results_df$CN == "Atypical"),
  Atypicality = which(results_df$Atypicality == "Atypical"),
  Norm = which(results_df$Norm == "Atypical")
)

# Find indices of samples identified as 'Atypical' by ALL methods
# Using Reduce with intersect to find the common intersection across the list
strong_outlier_indices <- Reduce(intersect, outlier_lists)

# Create a dataframe with only these extreme samples
strong_outliers_data <- results_df[strong_outlier_indices, ]

cat("--- EXTREME OUTLIER ANALYSIS ---\n")
cat("Number of samples identified by total consensus:", length(strong_outlier_indices), "\n")
cat("Percentage of total:", round(length(strong_outlier_indices)/N*100, 2), "%\n")

# Display the first indices
print(strong_outlier_indices)




### Investigating the nature of these outliers
# 1. Identify the "Strong" outliers
strong_idx <- Reduce(intersect, outlier_lists)

# 2. Create a factor in the original dataset for comparison
# Note: uses original 'Kola' dataset (non-ilr transformed) for chemical interpretation
Kola_analysis <- Kola
Kola_analysis$Type <- "Typical"
Kola_analysis$Type[strong_idx] <- "Strong Outlier"

# 3. Calculate geometric means (or medians) per group
# In geochemistry, the median is more robust for comparisons
library(dplyr)

# Geometric mean function (handles zeros with a small epsilon if necessary)
geom_mean <- function(x) {
  exp(mean(log(x[x > 0]), na.rm = TRUE))
}

# Calculation of the comparative table
comp_table_refined <- Kola_analysis %>%
  group_by(Type) %>%
  summarise(across(where(is.numeric), 
                   list(Med = median, GMean = geom_mean), 
                   .names = "{.col}_{.fn}")) %>%
  tidyr::pivot_longer(cols = -Type, 
                      names_to = c("Element", "Stat"), 
                      names_sep = "_") %>%
  tidyr::pivot_wider(names_from = Type, values_from = value) %>%
  mutate(Ratio = `Strong Outlier` / Typical) %>%
  arrange(Stat, desc(Ratio))

# Display Geometric Mean results for top elements
print("Most enriched elements (Geometric Mean):")
comp_table_refined %>% 
  filter(Stat == "GMean") %>% 
  head(10) %>% 
  print()

print("Top enriched elements in consensus outliers:")
print(head(comp_table_refined, 10))