## Does the unit-sphere constraint preserve the direction of beta when the
## data were generated with ||beta|| != 1?
##
## The identified quantities on the data side are the products
## (lam_q beta, lam_y beta, xi_1 beta) together with the loadings on the
## unit-variance latent noise. Restricting beta to the sphere therefore pins
## the scale, and the fitted beta should recover beta_true / ||beta_true||
## componentwise, with the loadings absorbing the norm. This script checks
## that at ||beta_true|| = 2 and 0.5.
##
## Run from the repository root:  Rscript simulation/test_beta_direction.R

source("R/lvc_fast.R")

M <- as.integer(Sys.getenv("LVC_M", "40000"))
gs <- 0.1

simfit <- function(cnorm) {
  set.seed(20260907)
  n <- 2000; p <- 4
  u_true <- c(0.6, -0.4, 0.3, 0.5); u_true <- u_true / sqrt(sum(u_true^2))
  beta_true <- cnorm * u_true
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

  fit <- run_lvc_fast(q, y, X, family = "Gaussian", M = M, burn = floor(M / 2),
                      gscale = gs, seed = 7, verbose = 0,
                      track_criteria = FALSE, beta_norm1 = TRUE)
  S <- fit$PAR[, (fit$burn + 1):fit$M]
  bhat <- rowMeans(S[paste0("beta", 1:p), ])
  blo  <- apply(S[paste0("beta", 1:p), ], 1, quantile, 0.025)
  bhi  <- apply(S[paste0("beta", 1:p), ], 1, quantile, 0.975)

  cat(sprintf("\n===== ||beta_true|| = %.2f =====\n", cnorm))
  tab <- cbind(u_true = u_true, post_mean = bhat, `2.5%` = blo, `97.5%` = bhi)
  rownames(tab) <- paste0("beta", 1:p)
  print(round(tab, 3))
  cat(sprintf("cosine(bhat, u_true) = %.5f\n", sum(bhat * u_true) /
                sqrt(sum(bhat^2) * sum(u_true^2))))
  cat("componentwise ratio bhat / u_true:",
      paste(round(bhat / u_true, 3), collapse = "  "), "\n")

  ## the loadings absorb the norm: lam_hat ~ cnorm * lam_true (approximately,
  ## since the unit latent variance also rescales the noise loading)
  lam <- rowMeans(S[c("lamq", "lamy", "xi1"), ])
  cat(sprintf("lamq: %.3f (true %.3f, x%.2f)   lamy: %.3f (true %.3f, x%.2f)   xi1: %.2f (true %.2f, x%.2f)\n",
              lam[1], lamq, lam[1] / lamq, lam[2], lamy, lam[2] / lamy,
              lam[3], xi[2], lam[3] / xi[2]))
  invisible(NULL)
}

simfit(2)
simfit(0.5)
