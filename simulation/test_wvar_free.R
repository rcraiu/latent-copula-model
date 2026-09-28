## Is the normalisation var(W | x) = 1 still needed once ||beta|| = 1?
##
## With beta on the unit sphere AND a free conditional variance s2W, the model
## is an exact reparametrisation of the free-beta model: the data identify
## lam * beta (hence lam and the direction) and lam * sd(W noise) (hence s2W).
## So the constrained fit should recover, for data generated with
## W ~ N(beta0' x, 1) and ||beta0|| = C:
##     beta  -> beta0 / C          (direction)
##     lamq  -> C * lamq0,  lamy -> C * lamy0,  xi1 -> C * xi10   (exactly)
##     s2W   -> 1 / C^2
## against the approximate multipliers seen when s2W is fixed at 1.
##
## Run from the repository root:  Rscript simulation/test_wvar_free.R

source("R/lvc_gamma_gauss.R")
suppressMessages(library(coda))

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

  fit <- run_lvc(q, y, X, family = "Gaussian", M = M, burn = floor(M / 2),
                 gscale = gs, seed = 7, verbose = 0,
                 track_criteria = FALSE, beta_norm1 = TRUE, wvar_free = TRUE)
  S <- fit$PAR[, (fit$burn + 1):fit$M]
  bhat <- rowMeans(S[paste0("beta", 1:p), ])
  blo  <- apply(S[paste0("beta", 1:p), ], 1, quantile, 0.025)
  bhi  <- apply(S[paste0("beta", 1:p), ], 1, quantile, 0.975)

  cat(sprintf("\n===== ||beta_true|| = %.2f, free var(W|x) =====\n", cnorm))
  tab <- cbind(u_true = u_true, post_mean = bhat, `2.5%` = blo, `97.5%` = bhi)
  rownames(tab) <- paste0("beta", 1:p)
  print(round(tab, 3))
  covered <- (blo <= u_true) & (u_true <= bhi)
  cat(sprintf("cosine(bhat, u_true) = %.5f   components covered: %d/4\n",
              sum(bhat * u_true) / sqrt(sum(bhat^2)), sum(covered)))

  lam <- rowMeans(S[c("lamq", "lamy", "xi1"), ])
  cat(sprintf("lamq: %.3f (expect %.3f)   lamy: %.3f (expect %.3f)   xi1: %.2f (expect %.2f)\n",
              lam[1], cnorm * lamq, lam[2], cnorm * lamy, lam[3], cnorm * xi[2]))
  s2 <- S["s2W", ]
  cat(sprintf("s2W: mean %.3f  95%% CI (%.3f, %.3f)   expect 1/C^2 = %.3f\n",
              mean(s2), quantile(s2, .025), quantile(s2, .975), 1 / cnorm^2))

  ## does freeing the variance open a new ridge?
  for (v in c("lamq", "lamy", "xi1"))
    cat(sprintf("cor(log s2W, %-4s) = %6.3f   ", v, cor(log(s2), S[v, ])))
  ess <- effectiveSize(mcmc(t(S[c("lamq", "lamy", "lalpha", "xi0", "xi1", "s2W",
                                  paste0("beta", 1:p)), seq(1, ncol(S), 5)])))
  cat("\nESS (thin 5):\n"); print(round(ess))
  cat(sprintf("acceptance: beta %.3f  W %.3f\n", fit$acc$beta, fit$acc$W))
  invisible(NULL)
}

simfit(2)
simfit(0.5)
simfit(1)
