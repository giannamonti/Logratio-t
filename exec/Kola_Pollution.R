# Reproducibility
set.seed(123)

source("exec/utils.R")       # confidence_ellipse_GM()
source("exec/fit_mvt_clr.R") # fit_mvt_clr(), dmt_clr(), maha_clr(),
# pseudo_log_det(), eff_rank()

## ---- Load packages ----
library(robCompositions)
library(compositions)
library(MASS)
library(mvtnorm)
library(fitHeavyTail)
library(mnormt)
library(StatDA)
library(Rfast)
library(ConfidenceEllipse)
library(caret)

## ---- Load data ----
data(chorizon)
Kola <- chorizon[, c("Co", "Cu", "Ni", "Mg", "Na", "S",
                     "As", "Bi", "Cd", "Sb", "Ag", "Pb")]
dim(Kola)

# =============================================================================
# CLR, ILR and ALR yield numerically identical Mahalanobis distances
# (max diff < 5e-14) and equivalent AIC values for distribution fitting.
# The model is formally defined in ILR coordinates (isometric, subcompositionally
# coherent). CLR is handled via ginv() and pseudo_log_det() in fit_mvt_clr().
# For 2-D visualisation of CLR (Panel 4, Figure 3), the fitted CLR mean and
# scatter are projected onto the 2-D SVD subspace — purely for plotting.
# =============================================================================

# SVD projection of CLR onto its (D-1) non-singular directions (plotting only)
clr_svd_project <- function(X_clr) {
  X_clr <- as.matrix(X_clr)
  sv    <- svd(scale(X_clr, scale = FALSE))
  V     <- sv$v[, seq_len(ncol(X_clr) - 1)]
  list(coords = X_clr %*% V, V = V)
}

###############################################################
# SECTION 1. Pollution subset  (D = 3, p = 2)
###############################################################

Pollution <- Kola[, c("Co", "Cu", "Ni")]
p <- 2L

# ---- ILR --------------------------------------------------------------------
Xcoord    <- as.matrix(ilr(Pollution))

fit_Norm  <- mvnorm.mle(Xcoord)
mu1       <- fit_Norm$mu;  sigma1 <- fit_Norm$sigma
loglik    <- dmvnorm(Xcoord, mu1, sigma1, log = TRUE)
AIC.normal <- 2 * p - 2 * sum(loglik)

MLEM     <- fit_mvt(Xcoord, nu = "iterative", nu_iterative_method = "ECM")
mu       <- MLEM$mu;  gdl <- MLEM$nu;  sigma <- MLEM$scatter
hat.cov  <- sigma * gdl / (gdl - 2)
loglik.t <- dmt(Xcoord, mean = mu, S = sigma, df = gdl, log = TRUE)
AIC.t    <- 2 * (p + 1) - 2 * sum(loglik.t)

Xcoord_df    <- as.data.frame(Xcoord); colnames(Xcoord_df) <- c("x", "y")
ellipse_99   <- confidence_ellipse(Xcoord_df, x = x, y = y, conf_level = 0.99)
ellipse_95   <- confidence_ellipse(Xcoord_df, x = x, y = y, conf_level = 0.95)
ellipse_50   <- confidence_ellipse(Xcoord_df, x = x, y = y, conf_level = 0.50)
ellipse_99_t <- confidence_ellipse_GM(Xcoord_df, x = x, y = y, conf_level = 0.99, gdl = gdl)
ellipse_95_t <- confidence_ellipse_GM(Xcoord_df, x = x, y = y, conf_level = 0.95, gdl = gdl)
ellipse_50_t <- confidence_ellipse_GM(Xcoord_df, x = x, y = y, conf_level = 0.50, gdl = gdl)

cat("=== Section 1 — ILR ===\n")
cat("AIC Normal:", AIC.normal, "\nAIC Student-t:", AIC.t, "\nnu:", gdl, "\n")

# ---- ALR --------------------------------------------------------------------
Xcoord.alr <- as.matrix(alr(Pollution))

fit_Norm.alr   <- mvnorm.mle(Xcoord.alr)
mu1.alr        <- fit_Norm.alr$mu;  sigma1.alr <- fit_Norm.alr$sigma
loglik.alr     <- dmvnorm(Xcoord.alr, mu1.alr, sigma1.alr, log = TRUE)
AIC.normal.alr <- 2 * p - 2 * sum(loglik.alr)

MLEM.alr     <- fit_mvt(Xcoord.alr, nu = "iterative", nu_iterative_method = "ECM")
mu.alr       <- MLEM.alr$mu;  gdl.alr <- MLEM.alr$nu;  sigma.alr <- MLEM.alr$scatter
hat.cov.alr  <- sigma.alr * gdl.alr / (gdl.alr - 2)
loglik.t.alr <- dmt(Xcoord.alr, mean = mu.alr, S = sigma.alr, df = gdl.alr, log = TRUE)
AIC.t.alr    <- 2 * (p + 1) - 2 * sum(loglik.t.alr)

cat("\n=== Section 1 — ALR ===\n")
cat("AIC Normal:", AIC.normal.alr, "\nAIC Student-t:", AIC.t.alr, "\nnu:", gdl.alr, "\n")

# ---- CLR --------------------------------------------------------------------
# fit_mvt_clr() operates directly on the N x 3 singular CLR matrix.
# ginv() and pseudo_log_det() handle the rank-deficient covariance.
Xcoord.clr_raw <- as.matrix(clr(Pollution))

## Normal
mu1.clr      <- colMeans(Xcoord.clr_raw)
sigma1.clr   <- cov(Xcoord.clr_raw)
Sinv1.clr    <- ginv(sigma1.clr)
loglik.clr   <- -0.5 * (
  nrow(Xcoord.clr_raw) * (p * log(2 * pi) + pseudo_log_det(sigma1.clr)) +
    sum(maha_clr(Xcoord.clr_raw, mu1.clr, Sinv1.clr))
)
AIC.normal.clr <- 2 * p - 2 * loglik.clr

## Student-t
MLEM.clr     <- fit_mvt_clr(Xcoord.clr_raw)
mu.clr       <- MLEM.clr$mu;  gdl.clr <- MLEM.clr$nu;  sigma.clr <- MLEM.clr$scatter
hat.cov.clr  <- sigma.clr * gdl.clr / (gdl.clr - 2)
loglik.t.clr <- dmt_clr(Xcoord.clr_raw, mu.clr, sigma.clr, gdl.clr, log = TRUE)
AIC.t.clr    <- 2 * (p + 1) - 2 * sum(loglik.t.clr)

cat("\n=== Section 1 — CLR ===\n")
cat("AIC Normal:", AIC.normal.clr, "\nAIC Student-t:", AIC.t.clr, "\nnu:", gdl.clr, "\n")

## Summary table Section 1
aic_table_s1 <- data.frame(
  Transform  = c("CLR", "ILR", "ALR"),
  AIC_Normal = c(AIC.normal.clr, AIC.normal,     AIC.normal.alr),
  AIC_t      = c(AIC.t.clr,      AIC.t,          AIC.t.alr),
  nu         = c(gdl.clr,        gdl,            gdl.alr)
)
cat("\n=== AIC Table (Section 1 — Pollution) ===\n")
print(aic_table_s1)

# ---- Rotation check: CLR-SVD is a rigid rotation of ILR --------------------
# The SVD basis of CLR and the ILR basis are related by an orthogonal matrix R
# with det(R) = +1 (pure rotation). This confirms that the level curves in
# Panel 4 are a rotated version of those in Panel 2.
proj          <- clr_svd_project(Xcoord.clr_raw)
Xcoord.clr2d  <- proj$coords
mu.clr2d      <- as.numeric(t(proj$V) %*% mu.clr)
sigma.clr2d   <- t(proj$V) %*% sigma.clr %*% proj$V
sigma1.clr2d  <- t(proj$V) %*% sigma1.clr %*% proj$V

M  <- t(Xcoord) %*% Xcoord.clr2d
sv <- svd(M)
D_rot <- diag(c(1, det(sv$u %*% t(sv$v))))
R  <- sv$u %*% D_rot %*% t(sv$v)

cat("\n=== ILR -> CLR rotation ===\n")
cat("R'R (should be I2):\n");    print(round(t(R) %*% R, 8))
cat("det(R) (should be +1):\n"); print(round(det(R), 8))
cat("Rotation angle (degrees):", round(acos(R[1, 1]) * 180 / pi, 2), "\n")

Xcoord_clr_df <- as.data.frame(Xcoord.clr2d); colnames(Xcoord_clr_df) <- c("x", "y")
ellipse_99_clr   <- confidence_ellipse(Xcoord_clr_df, x = x, y = y, conf_level = 0.99)
ellipse_95_clr   <- confidence_ellipse(Xcoord_clr_df, x = x, y = y, conf_level = 0.95)
ellipse_50_clr   <- confidence_ellipse(Xcoord_clr_df, x = x, y = y, conf_level = 0.50)
ellipse_99_t_clr <- confidence_ellipse_GM(Xcoord_clr_df, x = x, y = y, conf_level = 0.99, gdl = gdl.clr)
ellipse_95_t_clr <- confidence_ellipse_GM(Xcoord_clr_df, x = x, y = y, conf_level = 0.95, gdl = gdl.clr)
ellipse_50_t_clr <- confidence_ellipse_GM(Xcoord_clr_df, x = x, y = y, conf_level = 0.50, gdl = gdl.clr)

## ---- Figure 3 ---------------------------------------------------------------
# Colour scheme (no in-panel legends — see figure caption):
#   red  (col=2)    — Normal  (ILR and CLR)
#   blue (col=4)    — Student-t ILR
#   seagreen4       — Student-t CLR

col_norm  <- 2
col_t_ilr <- 4
col_t_clr <- "seagreen4"

# pdf("Figure3_Kola.pdf", height = 8, width = 8)
par(mfrow = c(2, 2),
    mar   = c(2.5, 2.5, 1.8, 0.5),
    mgp   = c(1.4, 0.4, 0),
    tcl   = -0.25,
    oma   = c(0, 0, 0, 0))

# Panel 1 (top-left): Ternary — ILR
plot(acomp(Pollution), pch = 4, cex = 0.4, col = "grey60")
for (prob in c(0.5, 0.95, 0.99)) {
  r <- sqrt(qchisq(prob, df = p))
  ellipses(acomp(ilrInv(mu1)), ilrvar2clr(sigma1), r, col = col_norm,  lwd = 1.8)
}
for (prob in c(0.5, 0.95, 0.99)) {
  r <- sqrt(qf(prob, p, gdl) * p)
  ellipses(acomp(ilrInv(mu)),  ilrvar2clr(sigma),  r, col = col_t_ilr, lwd = 1.8)
}

# Panel 2 (top-right): ILR coordinate space
plot(Xcoord_df, cex = 0.4, col = "grey60", pch = 19,
     xlim = c(-1, 2), ylim = c(-1, 2),
     xlab = expression(h(x)[1]), ylab = expression(h(x)[2]))
lines(ellipse_50$x,      ellipse_50$y,      col = col_norm,  lwd = 1.8)
lines(ellipse_95$x,      ellipse_95$y,      col = col_norm,  lwd = 1.8)
lines(ellipse_99$x,      ellipse_99$y,      col = col_norm,  lwd = 1.8)
lines(ellipse_50_t[, 1], ellipse_50_t[, 2], col = col_t_ilr, lwd = 1.8)
lines(ellipse_95_t[, 1], ellipse_95_t[, 2], col = col_t_ilr, lwd = 1.8)
lines(ellipse_99_t[, 1], ellipse_99_t[, 2], col = col_t_ilr, lwd = 1.8)

# Panel 3 (bottom-left): Ternary — CLR
plot(acomp(Pollution), pch = 4, cex = 0.4, col = "grey60")
for (prob in c(0.5, 0.95, 0.99)) {
  r <- sqrt(qchisq(prob, df = p))
  ellipses(acomp(clrInv(mu1.clr)), sigma1.clr, r, col = col_norm,  lwd = 1.8)
}
for (prob in c(0.5, 0.95, 0.99)) {
  r <- sqrt(qf(prob, p, gdl.clr) * p)
  ellipses(acomp(clrInv(mu.clr)), sigma.clr,   r, col = col_t_clr, lwd = 1.8)
}

# Panel 4 (bottom-right): CLR 2-D SVD projection (visualisation only)
all_x <- c(Xcoord_clr_df$x, ellipse_99_clr$x,      ellipse_99_t_clr[, 1])
all_y <- c(Xcoord_clr_df$y, ellipse_99_clr$y,      ellipse_99_t_clr[, 2])
xlim  <- range(all_x, na.rm = TRUE) * 1.05
ylim  <- range(all_y, na.rm = TRUE) * 1.05

plot(Xcoord_clr_df, cex = 0.4, col = "grey60", pch = 19,
     xlim = xlim, ylim = ylim, asp = 1,
     xlab = "CLR-PC1", ylab = "CLR-PC2")
lines(ellipse_50_clr$x,      ellipse_50_clr$y,      col = col_norm,  lwd = 1.8)
lines(ellipse_95_clr$x,      ellipse_95_clr$y,      col = col_norm,  lwd = 1.8)
lines(ellipse_99_clr$x,      ellipse_99_clr$y,      col = col_norm,  lwd = 1.8)
lines(ellipse_50_t_clr[, 1], ellipse_50_t_clr[, 2], col = col_t_clr, lwd = 1.8)
lines(ellipse_95_t_clr[, 1], ellipse_95_t_clr[, 2], col = col_t_clr, lwd = 1.8)
lines(ellipse_99_t_clr[, 1], ellipse_99_t_clr[, 2], col = col_t_clr, lwd = 1.8)
# dev.off()

###############################################################
# SECTION 2. Whole composition  (D = 12, p = 11)
###############################################################

p2 <- 11L

## ---- ILR -------------------------------------------------------------------
Xcoord2      <- as.matrix(ilr(Kola))

fit_Norm2    <- mvnorm.mle(Xcoord2)
mu1.2        <- fit_Norm2$mu;  sigma1.2 <- fit_Norm2$sigma
loglik2      <- dmvnorm(Xcoord2, mu1.2, sigma1.2, log = TRUE)
AIC.normal2  <- 2 * p2 - 2 * sum(loglik2)

MLEM2        <- fit_mvt(Xcoord2, nu = "iterative", nu_iterative_method = "ECM")
mu2          <- MLEM2$mu;  gdl2 <- MLEM2$nu;  sigma2 <- MLEM2$scatter
loglik.t2    <- dmt(Xcoord2, mean = mu2, S = sigma2, df = gdl2, log = TRUE)
AIC.t2       <- 2 * (p2 + 1) - 2 * sum(loglik.t2)

cat("\n=== Section 2 — ILR ===\n")
cat("AIC Normal:", AIC.normal2, "\nAIC Student-t:", AIC.t2, "\nnu:", gdl2, "\n")

## ---- ALR -------------------------------------------------------------------
Xcoord.alr2    <- as.matrix(na.omit(alr(Kola)))

fit_Norm.alr2  <- mvnorm.mle(Xcoord.alr2)
mu1.alr2       <- fit_Norm.alr2$mu;  sigma1.alr2 <- fit_Norm.alr2$sigma
loglik.alr2    <- dmvnorm(Xcoord.alr2, mu1.alr2, sigma1.alr2, log = TRUE)
AIC.normal.alr2 <- 2 * p2 - 2 * sum(loglik.alr2)

MLEM.alr2      <- fit_mvt(Xcoord.alr2, nu = "iterative", nu_iterative_method = "ECM")
mu.alr2        <- MLEM.alr2$mu;  gdl.alr2 <- MLEM.alr2$nu;  sigma.alr2 <- MLEM.alr2$scatter
loglik.t.alr2  <- dmt(Xcoord.alr2, mean = mu.alr2, S = sigma.alr2, df = gdl.alr2, log = TRUE)
AIC.t.alr2     <- 2 * (p2 + 1) - 2 * sum(na.omit(loglik.t.alr2))

cat("\n=== Section 2 — ALR ===\n")
cat("AIC Normal:", AIC.normal.alr2, "\nAIC Student-t:", AIC.t.alr2, "\nnu:", gdl.alr2, "\n")

## ---- CLR -------------------------------------------------------------------
Xcoord.clr2_raw <- as.matrix(clr(Kola))

## Normal
mu1.clr2      <- colMeans(Xcoord.clr2_raw)
sigma1.clr2   <- cov(Xcoord.clr2_raw)
Sinv1.clr2    <- ginv(sigma1.clr2)
loglik.clr2   <- -0.5 * (
  nrow(Xcoord.clr2_raw) * (p2 * log(2 * pi) + pseudo_log_det(sigma1.clr2)) +
    sum(maha_clr(Xcoord.clr2_raw, mu1.clr2, Sinv1.clr2))
)
AIC.normal.clr2 <- 2 * p2 - 2 * loglik.clr2

## Student-t
MLEM.clr2      <- fit_mvt_clr(Xcoord.clr2_raw)
mu.clr2        <- MLEM.clr2$mu;  gdl.clr2 <- MLEM.clr2$nu;  sigma.clr2 <- MLEM.clr2$scatter
loglik.t.clr2  <- dmt_clr(Xcoord.clr2_raw, mu.clr2, sigma.clr2, gdl.clr2, log = TRUE)
AIC.t.clr2     <- 2 * (p2 + 1) - 2 * sum(loglik.t.clr2)

cat("\n=== Section 2 — CLR ===\n")
cat("AIC Normal:", AIC.normal.clr2, "\nAIC Student-t:", AIC.t.clr2, "\nnu:", gdl.clr2, "\n")

## Summary table Section 2 (CLR first, then ILR, ALR)
aic_table_s2 <- data.frame(
  Transform  = c("CLR", "ILR", "ALR"),
  AIC_Normal = c(AIC.normal.clr2, AIC.normal2,      AIC.normal.alr2),
  AIC_t      = c(AIC.t.clr2,      AIC.t2,           AIC.t.alr2),
  nu         = c(gdl.clr2,        gdl2,             gdl.alr2)
)
cat("\n=== AIC Table (Section 2 — Full Kola) ===\n")
print(aic_table_s2)

## ---- PCA biplot (Figure 4) -------------------------------------------------
X_clr <- clr(Kola)
pca   <- prcomp(X_clr, scale = TRUE)
scores_pca <- as.data.frame(pca$x[, 1:2]); colnames(scores_pca) <- c("x", "y")

ellipse_99_pca <- confidence_ellipse(scores_pca, x = x, y = y, conf_level = 0.99)
ellipse_95_pca <- confidence_ellipse(scores_pca, x = x, y = y, conf_level = 0.95)
ellipse_50_pca <- confidence_ellipse(scores_pca, x = x, y = y, conf_level = 0.50)

MLEM_pca <- fit_mvt(scores_pca, nu = "iterative", nu_iterative_method = "ECM")
gdl_pca  <- MLEM_pca$nu

ellipse_99_t_pca <- confidence_ellipse_GM(scores_pca, x = x, y = y, conf_level = 0.99, gdl = gdl_pca)
ellipse_95_t_pca <- confidence_ellipse_GM(scores_pca, x = x, y = y, conf_level = 0.95, gdl = gdl_pca)
ellipse_50_t_pca <- confidence_ellipse_GM(scores_pca, x = x, y = y, conf_level = 0.50, gdl = gdl_pca)

# pdf("Figure4_Kola_biplot_clear.pdf", height = 7, width = 7)
scaling_factor <- 8.5; label_factor <- 10.5

all_x <- c(scores_pca$x,
           ellipse_99_pca$x,      ellipse_99_t_pca[, 1],
           pca$rotation[, 1] * label_factor)
all_y <- c(scores_pca$y,
           ellipse_99_pca$y,      ellipse_99_t_pca[, 2],
           pca$rotation[, 2] * label_factor)

pad   <- 0.5
xlim  <- range(all_x, na.rm = TRUE) + c(-pad, pad)
ylim  <- range(all_y, na.rm = TRUE) + c(-pad, pad)
rng   <- max(diff(xlim), diff(ylim)) / 2
mid_x <- mean(xlim); mid_y <- mean(ylim)
xlim  <- c(mid_x - rng, mid_x + rng)
ylim  <- c(mid_y - rng, mid_y + rng)

par(mar = c(4, 4, 1, 1))
plot(scores_pca, cex = 0.5, col = "grey60", pch = 19,
     xlim = xlim, ylim = ylim,
     xlab = "Dimension 1", ylab = "Dimension 2")
lines(ellipse_50_pca$x,      ellipse_50_pca$y,      col = 2, lwd = 2)
lines(ellipse_95_pca$x,      ellipse_95_pca$y,      col = 2, lwd = 2)
lines(ellipse_99_pca$x,      ellipse_99_pca$y,      col = 2, lwd = 2)
lines(ellipse_50_t_pca[, 1], ellipse_50_t_pca[, 2], col = 4, lwd = 2)
lines(ellipse_95_t_pca[, 1], ellipse_95_t_pca[, 2], col = 4, lwd = 2)
lines(ellipse_99_t_pca[, 1], ellipse_99_t_pca[, 2], col = 4, lwd = 2)

arrows(0, 0,
       pca$rotation[, 1] * scaling_factor,
       pca$rotation[, 2] * scaling_factor,
       col = "darkred", lwd = 1.5, length = 0.1)
text(pca$rotation[, 1] * label_factor,
     pca$rotation[, 2] * label_factor,
     labels = rownames(pca$rotation),
     col = "darkred", cex = 0.9, font = 2)
# dev.off()