# =============================================================================
# Envelope approach — goodness-of-fit for Normal and Student-t
# Full-sample Mahalanobis distances compared with theoretical quantiles
# and simulation-based confidence bands (B = 999 replications).
#
# Figure layout (2x2):
#   Top-left:     Pollution (D=3)  — Normal
#   Top-right:    Pollution (D=3)  — Student-t  (xlim capped at 80)
#   Bottom-left:  Kola full (D=12) — Normal
#   Bottom-right: Kola full (D=12) — Student-t  (xlim capped at 150)
#
# Both ILR and CLR yield identical distances (max diff < 5e-14).
# ILR coordinates are used here for computational convenience.
# =============================================================================

library(compositions)
library(fitHeavyTail)
library(mvtnorm)
library(mnormt)
library(StatDA)
library(MASS)

set.seed(123)
B <- 999

data(chorizon)
Kola <- chorizon[, c("Co", "Cu", "Ni", "Mg", "Na", "S",
                     "As", "Bi", "Cd", "Sb", "Ag", "Pb")]

# -----------------------------------------------------------------------------
# Helper: envelope for multivariate Normal
# -----------------------------------------------------------------------------
envelope_normal <- function(X, B = 999) {
  X   <- as.matrix(X)
  N   <- nrow(X)
  d   <- ncol(X)
  fit <- mvnorm.mle(X)
  mu  <- fit$mu
  S   <- fit$sigma
  
  d2_obs   <- sort(mahalanobis(X, mu, S))
  pp       <- ppoints(N)
  q_theory <- qchisq(pp, df = d)
  
  env_mat <- replicate(B, {
    X_sim <- mvtnorm::rmvnorm(N, mean = mu, sigma = S)
    sort(mahalanobis(X_sim, colMeans(X_sim), cov(X_sim)))
  })
  env_low <- apply(env_mat, 1, quantile, 0.025)
  env_upp <- apply(env_mat, 1, quantile, 0.975)
  
  list(q_theory = q_theory, d2_obs = d2_obs,
       env_low = env_low, env_upp = env_upp, d = d)
}

# -----------------------------------------------------------------------------
# Helper: envelope for multivariate Student-t
# -----------------------------------------------------------------------------
envelope_t <- function(X, B = 999) {
  X   <- as.matrix(X)
  N   <- nrow(X)
  d   <- ncol(X)
  fit <- fit_mvt(X, nu = "iterative", nu_iterative_method = "ECM")
  mu  <- fit$mu
  S   <- fit$scatter
  nu  <- fit$nu
  
  d2_obs   <- sort(mahalanobis(X, mu, S))
  pp       <- ppoints(N)
  q_theory <- d * qf(pp, df1 = d, df2 = nu)
  
  env_mat <- replicate(B, {
    X_sim  <- mvtnorm::rmvt(N, delta = mu, sigma = S, df = nu)
    fit_i  <- fit_mvt(X_sim, nu = nu)   # nu fixed to reduce variability
    sort(mahalanobis(X_sim, fit_i$mu, fit_i$scatter))
  })
  env_low <- apply(env_mat, 1, quantile, 0.025)
  env_upp <- apply(env_mat, 1, quantile, 0.975)
  
  list(q_theory = q_theory, d2_obs = d2_obs,
       env_low = env_low, env_upp = env_upp, nu = nu, d = d)
}

# -----------------------------------------------------------------------------
# Helper: draw one envelope panel (with optional xlim cap)
# -----------------------------------------------------------------------------
draw_envelope <- function(env, title = "",
                          col_pts = "grey40",
                          col_env = "grey70",
                          col_diag = 2,
                          xmax = NULL) {
  
  q_th   <- env$q_theory
  d2     <- env$d2_obs
  env_lo <- env$env_low
  env_up <- env$env_upp
  
  # cap x axis if requested
  if (!is.null(xmax)) {
    keep   <- q_th <= xmax
    q_th   <- q_th[keep]
    d2     <- d2[keep]
    env_lo <- env_lo[keep]
    env_up <- env_up[keep]
  }
  
  xlim <- c(0, max(q_th))
  ylim <- c(0, max(c(d2, env_up), na.rm = TRUE) * 1.02)
  
  plot(q_th, d2,
       pch = 19, cex = 0.4, col = col_pts,
       xlim = xlim, ylim = ylim,
       xlab = "Theoretical quantiles",
       ylab = expression(Mahalanobis ~ d^2),
       main = title)
  polygon(c(q_th, rev(q_th)),
          c(env_lo, rev(env_up)),
          col = adjustcolor(col_env, alpha.f = 0.35), border = NA)
  lines(q_th, env_lo, col = col_env, lty = 2, lwd = 1.2)
  lines(q_th, env_up, col = col_env, lty = 2, lwd = 1.2)
  abline(0, 1, col = col_diag, lwd = 1.8)
}

# =============================================================================
# Section 1: Pollution subset  (D = 3, ILR -> 2-D)
# =============================================================================
Pollution <- Kola[, c("Co", "Cu", "Ni")]
Xilr_poll <- as.matrix(ilr(Pollution))

message("Computing envelopes for Pollution (D=3)...")
env_norm_poll <- envelope_normal(Xilr_poll, B = B)
env_t_poll    <- envelope_t(Xilr_poll,    B = B)
message("  nu (Student-t, Pollution): ", round(env_t_poll$nu, 3))

# =============================================================================
# Section 2: Full Kola composition  (D = 12, ILR -> 11-D)
# =============================================================================
Xilr_kola <- as.matrix(ilr(Kola))

message("Computing envelopes for Kola full (D=12)...")
env_norm_kola <- envelope_normal(Xilr_kola, B = B)
env_t_kola    <- envelope_t(Xilr_kola,    B = B)
message("  nu (Student-t, Kola): ", round(env_t_kola$nu, 3))

# =============================================================================
# Figure: 2x2 envelope plot
# =============================================================================

# pdf("Figure_Envelope.pdf", height = 8, width = 8)
par(mfrow = c(2, 2),
    mar   = c(3.5, 3.5, 2, 0.5),
    mgp   = c(2, 0.5, 0),
    tcl   = -0.25,
    oma   = c(0, 0, 0, 0))

# Panel 1 (top-left): Pollution — Normal (full range)
draw_envelope(env_norm_poll,                          # <-- era env_t_poll
              title = expression("Pollution — Normal ("*chi^2*")"))

# Panel 2 (top-right): Pollution — Student-t (capped at 80)
draw_envelope(env_t_poll,
              title = bquote("Pollution" ~ "-" ~ "Student-t" ~ 
                               (hat(nu) == .(round(env_t_poll$nu, 1)))),
              xmax = 80)

# Panel 3 (bottom-left): Kola full — Normal (full range)
draw_envelope(env_norm_kola,
              title = expression("Kola — Normal ("*chi^2*")"))

# Panel 4 (bottom-right): Kola full — Student-t (capped at 150)
draw_envelope(env_t_kola,
              title = bquote("Kola" ~ "-" ~ "Student-t" ~ 
                               (hat(nu) == .(round(env_t_kola$nu, 1)))),
              xmax = 150)
# dev.off()