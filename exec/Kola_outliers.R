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

# --- 1. Inizializzazione ---
Xcoord <- ilr(Kola)
N <- nrow(Xcoord)
d <- ncol(Xcoord)
n_df <- N - 1

# Vettori per memorizzare le distanze LOO
d2_loo_norm <- numeric(N)
d2_loo_t    <- numeric(N)
atyp_pvals  <- numeric(N)


for (i in 1:N) {
  X_menoi <- as.matrix(Xcoord[-i, ])
  z <- Xcoord[i, , drop = FALSE]
  
  # A. Normale LOO (per Atypicality Index e Distanza Normale LOO)
  fit_norm_i <- mvnorm.mle(X_menoi)
  mu_n <- fit_norm_i$mu
  sig_n <- fit_norm_i$sigma
  
  # Distanza per Atypicality (formula qy)
  qy <- (1/(1 + 1/N)) * mahalanobis(z, mu_n, sig_n)
  atyp_pvals[i] <- pbeta(qy / (qy + n_df), d/2, (n_df - d + 1)/2)
  
  # Distanza per Conf. Matrix Normale LOO
  d2_loo_norm[i] <- mahalanobis(z, mu_n, sig_n)
  
  # B. Student-t LOO
  fit_t_i <- fit_mvt(X_menoi, nu = "iterative", nu_iterative_method = "ECM")
  d2_loo_t[i] <- mahalanobis(z, fit_t_i$mu, fit_t_i$scatter)
  
  if(i %% 100 == 0) cat("Punto", i, "di", N, "\n")
}

# --- 2. Metodi Robusti (Full Sample - Non serve LOO qui) ---
mcd_res    <- covMcd(Xcoord)
mcd_scores <- mcd_res$mah

comed_res    <- covComed(Xcoord, reweight = TRUE)
comed_scores <- mahalanobis(Xcoord, comed_res$center, comed_res$cov)

fit_cn <- ContaminatedMixt::CNmixt(Xcoord, G = 1, contamination = TRUE, model = "VVV", verbose = FALSE)
pred_cn <- factor(ifelse(as.numeric(fit_cn$models[[1]]$v) < 0.5, "Atypical", "Normal"))


# --- 3. Predizioni e Soglie ---

# Soglia T (usiamo quello del modello Full per coerenza)
final_t_fit <- fit_mvt(Xcoord, nu = "iterative", nu_iterative_method = "ECM")
gdl_t <- final_t_fit$nu
thresh_t <- d * qf(0.95, d, gdl_t)
t_scores <- mahalanobis(Xcoord, final_t_fit$mu, final_t_fit$scatter)


nu_loo_values <- numeric(N)
# Eseguiamo un ciclo veloce solo per estrarre nu (senza ricalcolare tutto)
message("Estrazione valori nu via LOO...")
for (i in 1:N) {
  X_menoi <- as.matrix(Xcoord[-i, ])
  # Usiamo un'estrazione veloce del parametro nu
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


# --- 4. Verifica e Consenso ---
print("Conteggio Outlier (Con T-LOO ricalcolato):")
print(colSums(preds[,-1] == "Atypical"))

# Outlier Forti (Unanimità tra T_LOO, MCD, CN, COMCoDa, Norm e Atypicality)
strong_idx <- which(preds$T_LOO == "Atypical" & 
                      preds$MCD == "Atypical" & 
                      preds$CN == "Atypical" &
                      preds$COMCoDa == "Atypical" &
                      preds$Atypicality == "Atypical" & 
                      preds$Norm_LOO == "Atypical")

cat("Outlier identificati all'unanimità (incluso T-LOO):", length(strong_idx), "\n")

# --- 2. Sintesi Risultati ---

results_df <- data.frame(
  ID = 1:N,
  MCD = preds$MCD,
  COMCoDa = preds$COMCoDa,
  T_Student = preds$T_LOO,
  CN = preds$CN,
  Atypicality = preds$Atypicality,
  Norm =  preds$Norm_LOO
)

print("Conteggio Outlier Identificati:")
print(colSums(results_df[, -1] == "Atypical", na.rm = TRUE))



outlier_lists <- list(
  MCD = which(results_df$MCD == "Atypical"),
  COMCoDa = which(results_df$COMCoDa == "Atypical"),
  T_Student = which(results_df$T_Student == "Atypical"),
  CN = which(results_df$CN == "Atypical"),
  Atypicality = which(results_df$Atypicality == "Atypical"),
  Norm = which(results_df$Norm == "Atypical")
)

# Trova gli indici dei campioni identificati come 'Atypical' da TUTTI i 5 metodi
# Usiamo Reduce con intersect per trovare l'intersezione comune a tutta la lista
strong_outlier_indices <- Reduce(intersect, outlier_lists)

# Creiamo un dataframe con solo questi campioni estremi
strong_outliers_data <- results_df[strong_outlier_indices, ]

cat("--- ANALISI OUTLIER ESTREMI ---\n")
cat("Numero di campioni identificati all'unanimità:", length(strong_outlier_indices), "\n")
cat("Percentuale sul totale:", round(length(strong_outlier_indices)/N*100, 2), "%\n")

# Visualizza i primi indici
print(strong_outlier_indices)




### cerchiamo di capire la natura di questi outliers
# 1. Identifichiamo i 42 outlier "Strong"
strong_idx <- Reduce(intersect, outlier_lists)

# 2. Creiamo un fattore nel dataset originale per il confronto
# Nota: usa il dataset originale 'Kola' (non trasformato ilr) per l'interpretazione chimica
Kola_analysis <- Kola
Kola_analysis$Type <- "Typical"
Kola_analysis$Type[strong_idx] <- "Strong Outlier"

# 3. Calcoliamo le medie geometriche (o mediane) per gruppo
# In geochimica la mediana è più robusta per i confronti
library(dplyr)

# Funzione per la media geometrica (gestisce eventuali zeri con un piccolo epsilon se necessario)
geom_mean <- function(x) {
  exp(mean(log(x[x > 0]), na.rm = TRUE))
}

# Calcolo della tabella comparativa
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

# Visualizziamo solo i risultati della Media Geometrica per gli elementi top
print("Elementi più arricchiti (Media Geometrica):")
comp_table_refined %>% 
  filter(Stat == "GMean") %>% 
  head(10) %>% 
  print()

print("Elementi più arricchiti negli outlier unanimi:")
print(head(comp_table_refined, 10))
