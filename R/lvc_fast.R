## Rcpp-backed driver for the gamma-Gaussian latent variable copula model.
## Initialisation and output match run_lvc() in lvc_gamma_gauss.R; the MCMC loop
## itself runs in C++ (lvc_sampler.cpp).

suppressMessages(library(Rcpp))
source("R/copula_densities.R")

## must match enum Family in lvc_sampler.cpp
FAMILY_CODE <- c(Gaussian = 0L, Clayton = 1L, Frank = 2L, Gumbel = 3L,
                 Independence = 4L, Joe = 5L, SurvClayton = 6L)

if (!exists(".lvc_compiled")) {
  sourceCpp("R/lvc_sampler.cpp")
  .lvc_compiled <- TRUE
}

run_lvc_fast <- function(q, y, X, family = "Gaussian", M = 20000, burn = NULL,
                         gscale = 0.1, seed = 1729, verbose = 2000,
                         prior = list(s2_mu = 100, s2_lam = 100, s2_xi = 100,
                                      a_alpha = 0.01, b_alpha = 0.01,
                                      a_sig = 0.01, b_sig = 0.01,
                                      s2_beta = 1000),
                         track_criteria = TRUE, thin_W = 0,
                         blockwise = FALSE, cross_block = FALSE) {
  stopifnot(family %in% names(FAMILY_CODE))
  set.seed(seed)
  n <- length(q); pd <- ncol(X)
  if (is.null(burn)) burn <- floor(M / 2)

  ## ---- initialisation, identical to the R sampler ----------------------
  beta <- as.numeric(solve(crossprod(X), crossprod(X, y)))
  beta <- beta / sqrt(sum(beta^2))
  W    <- as.numeric(X %*% beta)

  fit_q <- glm(q ~ W, family = Gamma(link = "log"))
  muq0  <- coef(fit_q)[1]; lamq0 <- coef(fit_q)[2]
  al0   <- 1 / max(summary(fit_q)$dispersion, 1e-3)
  fit_y <- lm(y ~ W)
  muy0  <- coef(fit_y)[1]; lamy0 <- coef(fit_y)[2]
  sig0  <- sqrt(mean(resid(fit_y)^2))
  if (lamy0 < 0) { beta <- -beta; W <- -W; lamq0 <- -lamq0; lamy0 <- -lamy0 }

  p <- c(muq = as.numeric(muq0), lamq = as.numeric(lamq0), lalpha = log(al0),
         muy = as.numeric(muy0), lamy = as.numeric(lamy0), lsigma = log(sig0),
         xi0 = 0, xi1 = 0)

  if (family != "Independence") {         # profile-MLE start for xi
    obj <- function(v) {
      pp <- p; pp["xi0"] <- v[1]; pp["xi1"] <- v[2]
      -sum(lvc_logdens(q, y, W, pp, family, gscale))
    }
    op <- try(optim(c(0, 0), obj, method = "Nelder-Mead",
                    control = list(maxit = 400)), silent = TRUE)
    if (!inherits(op, "try-error") && is.finite(op$value)) p[c("xi0", "xi1")] <- op$par
  }

  bp_prec <- crossprod(X) + diag(1 / prior$s2_beta, pd)
  bp_cov  <- solve(bp_prec)
  bp_cholT <- t(chol(bp_cov))

  fit <- run_lvc_cpp(as.numeric(q), as.numeric(y), as.matrix(X),
                     FAMILY_CODE[[family]], as.integer(M), as.integer(burn),
                     gscale, as.numeric(p), as.numeric(beta),
                     bp_cov, bp_cholT, prior, as.integer(thin_W),
                     track_criteria, as.integer(verbose), blockwise, cross_block)

  rownames(fit$PAR) <- c(names(p), paste0("beta", seq_len(pd)))
  fit$family <- family
  if (thin_W <= 0) fit$Wsamp <- NULL
  if (!track_criteria) fit$criteria <- NULL
  fit
}
