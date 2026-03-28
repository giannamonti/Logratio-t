# if needed
# install.packages(c("compositions", "robustbase", "ICSOutlier", "ContaminatedMixt", "Surrogate", "pROC"))

library(compositions)     
library(robustbase)       # covMcd e covComed
library(ICSOutlier)       # ics.outlier
library(ContaminatedMixt) # CNmixt
library(Surrogate)        # function Randvec
library(pROC)             # to compute AUC
library(tictoc)
library(MASS)
library(mvtnorm)       
library(fitHeavyTail)
library(mnormt)
library(Rfast)

generate_comp_data <- function(n, p, alpha) {
  # 1. random mean
  mu_simplex <- as.numeric(Surrogate::RandVec(a=0, b=1, s=1, n=p, m=1)$RandVecOutput)
  mu_ilr <- ilr(mu_simplex)
  
  # 2. variance-covariance matrix (as in Divino et. al., 2026)
  if(p == 3) {
    Sigma_ilr <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  } else { 
    Sigma_ilr <- diag(rep(1, p-1))
    Sigma_ilr[Sigma_ilr == 0] <- 0.5
  }
  
  # 3. outlier proportion
  n_out <- floor(n * (1 - alpha))
  n_typ <- n - n_out
  
  # 4. ilr data (typical)
  typical_ilr <- MASS::mvrnorm(n_typ, mu_ilr, Sigma_ilr)
  
  # 5. outlier in ilr coordinates 
  outliers_ilr <- MASS::mvrnorm(n_out, mu_ilr, Sigma_ilr * 25)
  
  # Outliers that "fade" into the main distribution. 
  # Moderate shift + doubled variance for outliers.
  # outliers_ilr <- MASS::mvrnorm(n_out, mu_ilr + 4, Sigma_ilr * 4)
  
  # 6. Union and final transformation
  all_data_ilr <- rbind(typical_ilr, outliers_ilr)
  data_simplex <- ilrInv(all_data_ilr) # Converts the entire N x (p-1) block to N x p
  
  labels <- c(rep(0, n_typ), rep(1, n_out))
  
  return(list(x = data_simplex, y = labels))
}


compute_metrics <- function(real, pred, scores = NULL) {
  # Contingency Table: [Predicted, Actual]
  # We force c(0,1) levels to avoid errors if a method doesn't find outliers
  tab <- table(factor(pred, levels = c(0, 1)), factor(real, levels = c(0, 1)))
  
  TP <- tab[2, 2] # True Positive
  TN <- tab[1, 1] # True Negative
  FP <- tab[2, 1] # False Positive
  FN <- tab[1, 2] # False Negative
  
  sensitivity <- TP / (TP + FN) # or Recall
  specificity <- TN / (TN + FP)
  ppv <- ifelse((TP + FP) == 0, 0, TP / (TP + FP)) # Positive Predictive Value (Precision)
  npv <- ifelse((TN + FN) == 0, 0, TN / (TN + FN)) # Negative Predictive Value
  
  # AUC calculation (requires the pROC package)
  auc_val <- NA
  if (!is.null(scores)) {
    # We use try() to avoid crashes if the AUC calculation fails.
    auc_obj <- try(pROC::roc(real, scores, quiet = TRUE), silent = TRUE)
    if (!inherits(auc_obj, "try-error")) {
      auc_val <- as.numeric(pROC::auc(auc_obj))
    }
  }
  
  return(c(AUC = auc_val, Sens = sensitivity, Spec = specificity, PPV = ppv, NPV = npv))
}



compute_param_loo <- function(data) {
  # input: ilr-data
  
  n <- nrow(data)
  d <- ncol(data)
  
  d2_loo_norm  <- numeric(n)
  d2_loo_t     <- numeric(n)
  atyp_pvals   <- numeric(n)
  nu_estimates <- numeric(n)
  
  for (i in 1:n) {
    X_menoi <- as.matrix(data[-i, , drop = FALSE])
    z <- data[i, , drop = FALSE]
    
    # NORMAL LOO
    fit_n <- mvnorm.mle(X_menoi)
    mu_n  <- fit_n$mu
    sig_n <- fit_n$sigma
    
    d2_loo_norm[i] <- mahalanobis(z, mu_n, sig_n)
    
    # ATYPICALITY
    n_df <- n - 1
    qy <- (1 / (1 + 1/n)) * d2_loo_norm[i]
    atyp_pvals[i] <- pbeta(qy / (qy + n_df), d/2, (n_df - d + 1)/2)
    
    # STUDENT-T (ECM) LOO
    fit_t <- fit_mvt(X_menoi, nu = "iterative", nu_iterative_method = "ECM")
    d2_loo_t[i] <- mahalanobis(z, fit_t$mu, fit_t$scatter)
    nu_estimates[i] <- fit_t$nu
    
  }
  
  return(list(
    norm_scores = d2_loo_norm,
    atyp_pvals  = atyp_pvals,
    t_scores    = d2_loo_t,
    avg_nu      = mean(nu_estimates) # We return the average nu for the final threshold
  ))
}

# setting parameter simulation (as in Divino et al., 2026)

set.seed(123)
n_val <- c(100, 1000)
p_val <- c(3, 5)
alpha_val <- c(0.9, 0.6)
repetitions <- 100


tic()
results_list <- list()

all_sim_results <- data.frame()

for (n in n_val) {
  for (p in p_val) {
    for (alpha in alpha_val) {
      message(sprintf("Running: n=%d, p=%d, alpha=%.1f", n, p, alpha))
      
      for (r in 1:repetitions) {
        #  data generation
        sim <- generate_comp_data(n, p, alpha)
        data_ilr <- ilr(sim$x)
        
        # Simulation uses ILR coordinates throughout.
        # CLR and ALR yield identical Mahalanobis distances and outlier
        # classifications (max diff < 5e-14) — see equivalence checks in
        # outlier_analysis_clr.R. Running CLR here would double computation
        # time without adding information.
        
        # 1. MCD
        mcd_res <- covMcd(data_ilr)
        # We identify as outliers those with 0 weight or high robust distance
        # mcd_pred <- ifelse(mcd_res$mcd.wt == 0, 1, 0)
        mcd_scores <- mcd_res$mah # Mahalanobis distance
        # mcd_metrics <- compute_metrics(sim$y, mcd_pred, mcd_scores)
        mcd_pred_fixed <- ifelse(mcd_scores > qchisq(0.95, p-1), 1, 0) ## GM
        mcd_metrics <- compute_metrics(sim$y, mcd_pred_fixed, mcd_scores)
        
        # 2. COMCoDa
        comed_res <- covComed(data_ilr, reweight = TRUE)
        center_comed <- comed_res$center
        cov_comed <- comed_res$cov
        # comed_pred <- ifelse(as.numeric(comed_res$weights) == 0, 1, 0)
        comed_scores <- mahalanobis(data_ilr, center = center_comed, cov = cov_comed)
        # comed_metrics <- compute_metrics(sim$y, comed_pred, comed_scores)
        comed_pred_fixed <- ifelse(comed_scores > qchisq(0.95, p-1), 1, 0) ## GM
        comed_metrics <- compute_metrics(sim$y, comed_pred_fixed, comed_scores)
        
        # # 3. ICS
        # ics_metrics <- rep(NA, 5) # Default if it fails
        # ics_res <- tryCatch({
        #   fit <- ics2(data_ilr)
        #   out <- ics.outlier(fit)
        #   
        #   # If successful, calculate the metrics
        #   pred <- ifelse(out@outliers, 1, 0)
        #   scores <- out@ics.distances
        #   compute_metrics(sim$y, pred, scores)
        # }, error = function(e) return(rep(NA, 5)))
        # 
        # ics_metrics <- ics_res
        
        # 4. CN (Contaminated Normal)
        # Note: CNmixt requires the number of components (G)
        cn_metrics <- tryCatch({
          cn_res <- CNmixt(data_ilr, G = 1, contamination = TRUE, model = "VVV", 
                           parallel = FALSE, verbose = FALSE)
          best_mod <- cn_res$models[[1]]
          v_prob <- best_mod$v # probability of being 'good'
          if (is.null(v_prob)) {
            v_prob <- rep(1, nrow(data_ilr))
          } else {
            v_prob <- as.numeric(v_prob)
          }
          pred <- ifelse(v_prob < 0.5, 1, 0) # If the probability of being good is < 0.5, it is an outlier
          
          # The score for the AUC should be the probability of being an outlier
          scores <- 1 - v_prob
          
          compute_metrics(sim$y, pred, scores)
        }, error = function(e) {
          return(rep(NA, 5)) 
        })
        
        loo_results <- compute_param_loo(data_ilr)
        
        # thresholds and predictions
        d     <- ncol(data_ilr)
        n_obs <- nrow(data_ilr)  # renamed to avoid overwriting loop variable n
        
        norm_pred <- ifelse(loo_results$norm_scores > qchisq(0.95, d), 1, 0)
        atyp_pred <- ifelse(loo_results$atyp_pvals > 0.95, 1, 0)
        
        # Student-t prediction (using the mean of the estimated nu)
        thresh_t <- d * qf(0.95, d, loo_results$avg_nu)
        t_pred <- ifelse(loo_results$t_scores > thresh_t, 1, 0)
        
        norm_metrics <- compute_metrics(sim$y, norm_pred, loo_results$norm_scores)
        atyp_metrics <- compute_metrics(sim$y, atyp_pred, loo_results$atyp_pvals)
        t_metrics <- compute_metrics(sim$y, t_pred, loo_results$t_scores)
        
        temp_res <- rbind(t_metrics, norm_metrics, atyp_metrics,
                          mcd_metrics, comed_metrics, cn_metrics)
        rownames(temp_res) <- c("T","Norm", "Atyp","MCD", "COMCoDa", "CN")
        iteration_df <- as.data.frame(temp_res)
        iteration_df$Method <- rownames(temp_res)
        iteration_df$n <- n_obs
        iteration_df$p <- p
        iteration_df$alpha <- alpha
        iteration_df$rep <- r
        
        all_sim_results <- rbind(all_sim_results, iteration_df)
        print(r)
      }
    }
  }
}
toc()


library(dplyr)

final_table <- all_sim_results %>%
  group_by(n, p, alpha, Method) %>%
  summarise(
    Avg_AUC  = mean(AUC, na.rm = TRUE),
    Avg_Sens = mean(Sens, na.rm = TRUE),
    Avg_Spec = mean(Spec, na.rm = TRUE),
    Avg_PPV  = mean(PPV, na.rm = TRUE),
    Avg_NPV  = mean(NPV, na.rm = TRUE),
    .groups = "drop"
  )

final_table_formatted <- final_table %>%
  mutate(across(starts_with("Avg_"), ~ round(., 3)))

print(final_table_formatted)

# write.csv(final_table, "/Users/giannamonti/Downloads/risultati_simulazione_31dic25.csv")


library(ggplot2)
library(dplyr)

# 1. Data preparation: select the methods you want to compare
plot_data <- all_sim_results %>%
  filter(Method %in% c("Norm", "Atyp", "T", "T3", "MCD", "COMCoDa")) %>%
  mutate(alpha_label = paste0("Contamination: ", round((1-alpha)*100, 0), "%"))

# 2. Boxplot for Sensibility
ggplot(plot_data, aes(x = Method, y = Sens, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  facet_wrap(~ alpha_label) + 
  theme_minimal() +
  labs(
    title = " ",
    x = "Method",
    y = "Sensibility (True Positive Rate)",
    fill = "Method"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  )


