## Two PX-DA designs, both adapted from Xu, Craiu, Sun & Paterson (2016, eq. 3.2).
##
##   A. Paper's PX-HC: identify by var(W|x) = 1, beta free. Draw unrestricted
##      beta* and residual scale psi^2, reduce by psi.
##   B. Identified by ||beta|| = 1 with free s2W, plus a Metropolis group move
##      along the scale orbit (W, loadings, s2W).
##
## Run from the repository root:  Rscript simulation/test_px_da.R

if (exists(".lvc_compiled")) rm(.lvc_compiled)
suppressMessages(library(Rcpp))
sourceCpp("R/lvc_sampler.cpp")
.lvc_compiled <- TRUE
source("R/lvc_fast.R")
suppressMessages(library(coda))

M <- as.integer(Sys.getenv("LVC_M", "40000"))
gs <- 0.1

simdata <- function() {
  set.seed(20260907)
  n <- 2000; p <- 4
  u <- c(0.6, -0.4, 0.3, 0.5); u <- u / sqrt(sum(u^2))
  beta_true <- u                         # ||beta|| = 1, var(W|x) = 1
  muq <- 0; lamq <- 0.40; alpha <- 20
  muy <- -0.20; lamy <- 0.30; sig <- 0.25
  xi <- c(2, 4)
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X); X <- matrix(as.numeric(X), n, p)
  W <- as.numeric(X %*% beta_true) + rnorm(n)
  rho <- theta_of_w(W, xi, "Gaussian", gs)
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  q <- qgamma(pnorm(z1), shape = alpha, rate = alpha / exp(muq + lamq * W))
  y <- qnorm(pnorm(z2), mean = muy + lamy * W, sd = sig)
  list(q = q, y = y, X = X, u = u, lamq = lamq, lamy = lamy, xi1 = xi[2])
}

ess_of <- function(S, vars) {
  effectiveSize(mcmc(t(S[vars, seq(1, ncol(S), by = 5)])))
}

summarise <- function(fit, D, label, sphere = FALSE) {
  S <- fit$PAR[, (fit$burn + 1):fit$M]
  b <- S[paste0("beta", 1:4), ]
  bn <- sqrt(colSums(b^2))
  bhat <- rowMeans(b)
  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("cosine(bhat, u) = %.5f   ||beta|| mean %.3f  sd %.3f\n",
              sum(bhat * D$u) / sqrt(sum(bhat^2)), mean(bn), sd(bn)))
  cat(sprintf("lamq %.3f (true %.2f)  lamy %.3f (true %.2f)  xi1 %.2f (true %.2f)\n",
              mean(S["lamq", ]), D$lamq, mean(S["lamy", ]), D$lamy,
              mean(S["xi1", ]), D$xi1))
  if ("s2W" %in% rownames(S)) {
    s2 <- S["s2W", ]
    cat(sprintf("s2W mean %.3f  95%% CI (%.3f, %.3f)\n",
                mean(s2), quantile(s2, .025), quantile(s2, .975)))
    for (v in c("lamq", "lamy", "xi1"))
      cat(sprintf("cor(log s2W, %-4s) = %6.3f  ", v, cor(log(s2), S[v, ])))
    cat("\n")
  } else {
    cat(sprintf("cor(log||beta||, lamq) = %.3f   lamy = %.3f   xi1 = %.3f\n",
                cor(log(bn), S["lamq", ]), cor(log(bn), S["lamy", ]),
                cor(log(bn), S["xi1", ])))
  }
  vars <- c("lamq", "lamy", "lalpha", "xi0", "xi1", paste0("beta", 1:4))
  if ("s2W" %in% rownames(S)) vars <- c(vars, "s2W")
  ess <- ess_of(S, vars)
  cat("ESS (thin 5):\n"); print(round(ess))
  if (!is.null(fit$acc_px) && is.finite(fit$acc_px))
    cat(sprintf("group-move acceptance = %.3f\n", fit$acc_px))
  if (!is.null(fit$acc) && !is.null(fit$acc$px) && is.finite(fit$acc$px))
    cat(sprintf("group-move acceptance = %.3f\n", fit$acc$px))
  invisible(ess)
}

D <- simdata()

cat("---- A. identification var(W|x) = 1, beta free ----\n")
fitA0 <- run_lvc_fast(D$q, D$y, D$X, family = "Gaussian", M = M,
                      burn = floor(M / 2), gscale = gs, seed = 7, verbose = 0,
                      track_criteria = FALSE)
fitA1 <- run_lvc_fast(D$q, D$y, D$X, family = "Gaussian", M = M,
                      burn = floor(M / 2), gscale = gs, seed = 7, verbose = 0,
                      track_criteria = FALSE, px_da = TRUE)
eA0 <- summarise(fitA0, D, "unconstrained Gibbs (var = 1)")
eA1 <- summarise(fitA1, D, "PX-HC, eq. (3.2), reduce by psi")
cat(sprintf("\nESS ratio A (PX / unconstrained): lamq %.2f  lamy %.2f  xi1 %.2f  beta1 %.2f\n",
            eA1["lamq"] / eA0["lamq"], eA1["lamy"] / eA0["lamy"],
            eA1["xi1"] / eA0["xi1"], eA1["beta1"] / eA0["beta1"]))

cat("\n---- B. identification ||beta|| = 1, s2W free ----\n")
fitB0 <- run_lvc_fast(D$q, D$y, D$X, family = "Gaussian", M = M,
                      burn = floor(M / 2), gscale = gs, seed = 7, verbose = 0,
                      track_criteria = FALSE, beta_norm1 = TRUE, wvar_free = TRUE)
fitB1 <- run_lvc_fast(D$q, D$y, D$X, family = "Gaussian", M = M,
                      burn = floor(M / 2), gscale = gs, seed = 7, verbose = 0,
                      track_criteria = FALSE, beta_norm1 = TRUE, wvar_free = TRUE,
                      px_da = TRUE)
eB0 <- summarise(fitB0, D, "naive free-variance Gibbs")
eB1 <- summarise(fitB1, D, "free-variance + PX group move")
cat(sprintf("\nESS ratio B (group / naive): lamq %.2f  lamy %.2f  xi1 %.2f  s2W %.2f\n",
            eB1["lamq"] / eB0["lamq"], eB1["lamy"] / eB0["lamy"],
            eB1["xi1"] / eB0["xi1"], eB1["s2W"] / eB0["s2W"]))
