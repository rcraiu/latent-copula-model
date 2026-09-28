## Does constraining ||beta|| = 1 resolve the scale ridge?
##
## The likelihood is (nearly) invariant under W -> W/C, beta -> beta/C,
## (lam_q, lam_y, xi_1) -> C (lam_q, lam_y, xi_1): only the N(beta'x, 1)
## prior on W resists the rescaling, so the posterior has a flat ridge
## trading the norm of beta against the loading parameters. This script
## fits the recovery scenario of test_recovery.R (where short runs fail
## to recover all parameters) twice on the same data -- once with the
## unconstrained conjugate beta draw, once with beta restricted to the
## unit sphere -- and compares ridge evidence, mixing, and recovery.
##
## Run from the repository root:  Rscript simulation/test_beta_norm1.R

source("R/lvc_gamma_gauss.R")

ess <- function(v) {                     # coda if present, else batch means
  if (requireNamespace("coda", quietly = TRUE))
    return(as.numeric(coda::effectiveSize(v)))
  nb <- 50; m <- floor(length(v) / nb)
  bm <- colMeans(matrix(v[1:(nb * m)], m, nb))
  length(v) * var(v) / (m * var(bm) * length(v) / (m * nb)) # ~ n * s2 / spec0
}

set.seed(20260907)                       # identical data to test_recovery.R
n  <- 2000
p  <- 4
gs <- 0.1

beta_true <- c(0.6, -0.4, 0.3, 0.5); beta_true <- beta_true / sqrt(sum(beta_true^2))
muq_true  <- 0.00; lamq_true <- 0.40; alpha_true <- 20
muy_true  <- -0.20; lamy_true <- 0.30; sig_true  <- 0.25
xi_true   <- c(2, 4)

X <- matrix(rnorm(n * p), n, p)
X <- scale(X); X <- matrix(as.numeric(X), n, p)
W <- as.numeric(X %*% beta_true) + rnorm(n)

rho <- theta_of_w(W, xi_true, "Gaussian", gs)
z1  <- rnorm(n)
z2  <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
u1  <- pnorm(z1); u2 <- pnorm(z2)

mq <- exp(muq_true + lamq_true * W)
q  <- qgamma(u1, shape = alpha_true, rate = alpha_true / mq)
y  <- qnorm(u2, mean = muy_true + lamy_true * W, sd = sig_true)

M <- as.integer(Sys.getenv("LVC_M", "40000"))

truth <- c(muq = muq_true, lamq = lamq_true, lalpha = log(alpha_true),
           muy = muy_true, lamy = lamy_true, lsigma = log(sig_true),
           xi0 = xi_true[1], xi1 = xi_true[2],
           beta1 = beta_true[1], beta2 = beta_true[2],
           beta3 = beta_true[3], beta4 = beta_true[4])

report <- function(fit, label) {
  S  <- fit$PAR[, (fit$burn + 1):fit$M]
  bn <- sqrt(colSums(S[paste0("beta", 1:p), ]^2))
  cat("\n============", label, "============\n")

  ## ridge evidence: norm of beta against the loadings it trades against
  cat(sprintf("||beta||: mean %.3f  sd %.4f\n", mean(bn), sd(bn)))
  for (v in c("lamq", "lamy", "xi1"))
    cat(sprintf("  cor(log||beta||, %-4s) = %+.3f\n", v, cor(log(bn), S[v, ])))

  ## mixing and recovery for the scale-coupled parameters
  keep <- names(truth)
  tab <- t(sapply(keep, function(v) {
    dr <- S[v, ]
    c(truth = truth[[v]], mean = mean(dr),
      lo = quantile(dr, 0.025), hi = quantile(dr, 0.975), ESS = ess(dr))
  }))
  colnames(tab) <- c("truth", "mean", "2.5%", "97.5%", "ESS")
  print(round(tab, 3))
  cov <- tab[, "truth"] >= tab[, "2.5%"] & tab[, "truth"] <= tab[, "97.5%"]
  cat("covered:", sum(cov), "of", length(cov),
      if (any(!cov)) paste0("  (missed: ", paste(keep[!cov], collapse = ", "), ")") else "", "\n")
  cat("ESS(||beta||) =", round(ess(bn), 1), "\n")
  invisible(tab)
}

cat("M =", M, "iterations, burn =", floor(M / 2), "\n")

t1 <- system.time(
  fit_u <- run_lvc(q, y, X, family = "Gaussian", M = M, burn = floor(M / 2),
                   gscale = gs, seed = 7, verbose = 0, track_criteria = FALSE)
)
cat("\nunconstrained run:", round(t1[["elapsed"]], 1), "sec\n")
tab_u <- report(fit_u, "UNCONSTRAINED (conjugate beta)")

t2 <- system.time(
  fit_c <- run_lvc(q, y, X, family = "Gaussian", M = M, burn = floor(M / 2),
                   gscale = gs, seed = 7, verbose = 0, track_criteria = FALSE,
                   beta_norm1 = TRUE)
)
cat("\nconstrained run:", round(t2[["elapsed"]], 1), "sec\n")
tab_c <- report(fit_c, "CONSTRAINED (||beta|| = 1)")
cat("\nbeta step acceptance rate:", round(fit_c$acc$beta, 3), "\n")

cat("\nESS ratio (constrained / unconstrained):\n")
print(round(tab_c[, "ESS"] / tab_u[, "ESS"], 2))
