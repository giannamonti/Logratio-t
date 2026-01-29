# Reproducibility: set a fixed random seed

set.seed(123)

source("exec/utils.R")

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

###############################################################
# SECTION 1. Pollution subset
###############################################################

Pollution <- Kola[, c("Co", "Cu", "Ni")]

# ILR Analysis
V <- matrix(c(-1/sqrt(2), 1/sqrt(2), 0, 
              -1/sqrt(6), -1/sqrt(6), sqrt(2/3)), ncol = 2)
Xcoord <- ilr(Pollution)
Xcoord <- as.matrix(Xcoord, drop = FALSE)

## ---- Fit multivariate normal ----
fit_Norm <- mvnorm.mle(Xcoord)
mu1 <- fit_Norm$mu
sigma1 <- fit_Norm$sigma
n <- nrow(Xcoord); p <- length(mu1)
loglik <- dmvnorm(Xcoord, mu1, sigma1, log = TRUE)
AIC.normal <- 2 * p - 2 * sum(loglik)

## ---- Confidence ellipses (normal) ----
Xcoord_df <- as.data.frame(Xcoord)
colnames(Xcoord_df) <- c("x", "y")
ellipse_99 <- confidence_ellipse(Xcoord_df, x = x, y = y, conf_level = 0.99)
ellipse_95 <- confidence_ellipse(Xcoord_df, x = x, y = y, conf_level = 0.95)
ellipse_50 <- confidence_ellipse(Xcoord_df, x = x, y = y, conf_level = 0.50)

## ---- Fit multivariate Student-t ----
MLEM <- fit_mvt(Xcoord, nu = "iterative", nu_iterative_method = "ECM")
mu <- MLEM$mu; gdl <- MLEM$nu; sigma <- MLEM$scatter; hat.cov <- sigma*gdl/(gdl - 2)
loglik.t <- dmt(Xcoord, mean = mu, S = sigma, df = gdl, log = TRUE)
AIC.t <- 2 * (p + 1) - 2 * sum(loglik.t)

ellipse_99_t <- confidence_ellipse_GM(Xcoord_df, x = x, y = y, conf_level = 0.99, gdl = gdl)
ellipse_95_t <- confidence_ellipse_GM(Xcoord_df, x = x, y = y, conf_level = 0.95, gdl = gdl)
ellipse_50_t <- confidence_ellipse_GM(Xcoord_df, x = x, y = y, conf_level = 0.50, gdl = gdl)


# ALR Analysis

Xcoord.alr <- alr(Pollution)
Xcoord.alr <- as.matrix(Xcoord.alr, drop = FALSE)
MLEM.alr<- fit_mvt(Xcoord.alr, nu = "iterative", nu_iterative_method = "ECM")
mu.alr <- MLEM.alr$mu; gdl.alr <- MLEM.alr$nu; sigma.alr <- MLEM.alr$scatter; hat.cor.alr <- sigma.alr*gdl.alr/(gdl.alr - 2)
loglik.t.alr <- dmt(as.matrix(Xcoord.alr), mean=MLEM.alr$mu,S=sigma.alr,df=gdl.alr, log = T)
AIC.t.alr <- 2*(p+1) - 2*sum(loglik.t.alr)

fit_Norm.alr <- mvnorm.mle(Xcoord.alr)
mu1.alr <- fit_Norm.alr$mu
sigma1.alr <- fit_Norm.alr$sigma
loglik.alr <- dmvnorm(Xcoord.alr, mu1.alr, sigma1.alr, log = TRUE)
AIC.normal.alr <- 2 * p - 2 * sum(loglik.alr)


## ---- Figure 3: Pollution subset ----
# pdf("Figure3_Kola.pdf", height = 5, width = 9)
par(mfrow = c(1, 2))

plot(acomp(Pollution), pch = 4, cex = 0.5, col = "grey60")
for (p in c(0.5, 0.95, 0.99)) {
  r = sqrt(qchisq(p = p, df = 2))
  ellipses(acomp(ilrInv(mu1)), ilrvar2clr(sigma1), r, col = 2, lwd = 2)
}
for (p in c(0.5, 0.95, 0.99)) {
  r = sqrt(qf(p, 2, gdl) * 2)
  ellipses(acomp(ilrInv(mu)), ilrvar2clr(sigma), r, col = 4, lwd = 2)
}

plot(Xcoord_df, cex = .5, col = "grey60", pch = 19,
     xlim = c(-1, 2), ylim = c(-1, 2), 
     xlab = expression(h(x)[1]), ylab = expression(h(x)[2]))
lines(ellipse_50$x, ellipse_50$y, col = 2, lwd = 2)
lines(ellipse_95$x, ellipse_95$y, col = 2, lwd = 2)
lines(ellipse_99$x, ellipse_99$y, col = 2, lwd = 2)
lines(ellipse_50_t[, 1], ellipse_50_t[, 2], col = 4, lwd = 2)
lines(ellipse_95_t[, 1], ellipse_95_t[, 2], col = 4, lwd = 2)
lines(ellipse_99_t[, 1], ellipse_99_t[, 2], col = 4, lwd = 2)
# dev.off()

###############################################################
# SECTION 2. Whole composition
###############################################################

## ---- Likelihood-based goodness of fit measures (Table 2) ----

Xcoord <- ilr(Kola)
fit_Norm <- mvnorm.mle(as.matrix(Xcoord))
mu1 <- fit_Norm$mu
sigma1 <- fit_Norm$sigma
loglik <- dmvnorm(as.matrix(Xcoord), mu1, sigma1, log = TRUE)
p <- length(mu1); n <- nrow(Xcoord)
AIC.normal <- 2 * p - 2 * sum(loglik)

MLEM <- fit_mvt(Xcoord, nu = "iterative", nu_iterative_method = "ECM")
mu <- MLEM$mu; gdl <- MLEM$nu; sigma <- MLEM$scatter
loglik.t <- dmt(Xcoord, mean = mu, S = sigma, df = gdl, log = TRUE)
AIC.t <- 2 * (p + 1) - 2 * sum(loglik.t)


Xcoord.alr <- alr(Kola)
Xcoord.alr <-na.omit(Xcoord.alr)
MLEM.alr <- fit_mvt(Xcoord.alr, nu = "iterative", nu_iterative_method = "ECM")
mu.alr <- MLEM.alr$mu; gdl.alr <- MLEM.alr$nu; sigma.alr <- MLEM.alr$scatter
loglik.t.alr <- dmt(as.matrix(Xcoord.alr), mean=mu.alr,S=sigma.alr,df=gdl.alr, log = T)
AIC.t.alr <- 2*(p+1) - 2*sum(na.omit(loglik.t.alr))

fit_Norm.alr <- mvnorm.mle(as.matrix(Xcoord.alr))
mu1.alr <- fit_Norm.alr$mu
sigma1.alr <- fit_Norm.alr$sigma
loglik.alr <- dmvnorm(as.matrix(Xcoord.alr), mu1.alr, sigma1.alr, log = TRUE)
AIC.normal.alr <- 2 * p - 2 * sum(loglik.alr)


## ---- PCA and biplot (Figure 4) ----
X <- clr(Kola)
pca <- prcomp(X, scale = TRUE)
scores <- pca$x[, 1:2]
Xcoord <- as.data.frame(scores)
colnames(Xcoord) <- c("x", "y")

ellipse_99 <- confidence_ellipse(Xcoord, x = x, y = y, conf_level = 0.99)
ellipse_95 <- confidence_ellipse(Xcoord, x = x, y = y, conf_level = 0.95)
ellipse_50 <- confidence_ellipse(Xcoord, x = x, y = y, conf_level = 0.50)
MLEM <- fit_mvt(Xcoord, nu = "iterative", nu_iterative_method = "ECM")
gdl <- MLEM$nu
ellipse_99_t <- confidence_ellipse_GM(Xcoord, x = x, y = y, conf_level = 0.99, gdl = gdl)
ellipse_95_t <- confidence_ellipse_GM(Xcoord, x = x, y = y, conf_level = 0.95, gdl = gdl)
ellipse_50_t <- confidence_ellipse_GM(Xcoord, x = x, y = y, conf_level = 0.50, gdl = gdl)


# pdf("Figure4_Kola_biplot_clear.pdf", height = 7, width = 7)
par(pty = "s", mar = c(5, 5, 4, 4))

plot(Xcoord, cex = .5, col = "grey60", pch = 19,
     xlim = c(-12, 12), ylim = c(-12, 12), 
     xlab = "Dimension 1", ylab = "Dimension 2")

lines(ellipse_50$x, ellipse_50$y, col = 2, lwd = 2)
lines(ellipse_95$x, ellipse_95$y, col = 2, lwd = 2)
lines(ellipse_99$x, ellipse_99$y, col = 2, lwd = 2)
lines(ellipse_50_t[, 1], ellipse_50_t[, 2], col = 4, lwd = 2)
lines(ellipse_95_t[, 1], ellipse_95_t[, 2], col = 4, lwd = 2)
lines(ellipse_99_t[, 1], ellipse_99_t[, 2], col = 4, lwd = 2)

scaling_factor <- 8.5
arrows(0, 0, pca$rotation[, 1] * scaling_factor, pca$rotation[, 2] * scaling_factor, 
       col = "darkred", lwd = 1.5, length = 0.1)
label_factor <- 10.5 
text(pca$rotation[, 1] * label_factor, 
     pca$rotation[, 2] * label_factor, 
     labels = rownames(pca$rotation),
     col = "darkred", 
     cex = 0.9, 
     font = 2)
# dev.off()


