## Model functions for the mixed-response latent variable copula model
## (Poisson count margin, normal margin), Section 2.2 of the paper:
##
##   Q_i | W_i ~ Poisson( exp(mu_d + lam_d W_i) )
##   Y_i | W_i ~ N( mu_c + lam_c W_i, sigma^2 )
##   (Q_i, Y_i) | W_i ~ C( F_d, F_c | theta(W_i, xi) ),  theta = kappa(gam (xi_0 + xi_1 W_i))
##   W_i ~ N( beta' x_i, 1 )
##
## Everything is vectorised over observations with the copula parameter varying
## across them. The activation functions and copula log densities are shared
## with the continuous-continuous code (R/copula_densities.R). Run from the
## repository root.
##
## Notation: h(u1 | u2; theta) = dC(u1, u2)/du2 = P(U1 <= u1 | U2 = u2) is the
## conditional copula C_{d|c} of the paper. The observed-data likelihood is
##   p(q, y | w, eta) = [ h(F_d(q) | F_c(y)) - h(F_d(q - 1) | F_c(y)) ] f_c(y | w).

source("R/copula_densities.R")

## ---------------------------------------------------------------- h-functions
## The parameter is clamped to a range in which the closed forms are finite in
## double precision; the ranges are far wider than the posterior ever visits.
clamp_theta <- function(th, family) {
  switch(family,
    Gaussian     = pmax(pmin(th, 1 - 1e-8), -1 + 1e-8),
    Clayton      = pmax(pmin(th, 200), 1e-8),
    Gumbel       = pmax(pmin(th, 200), 1 + 1e-8),
    Frank        = { th <- pmax(pmin(th, 200), -200); ifelse(abs(th) < 1e-8, 1e-8, th) },
    Independence = th,
    stop("unknown family: ", family))
}

hfunc <- function(u1, u2, th, family) {
  u2 <- clamp01(u2); v <- clamp01(u1)
  th <- clamp_theta(th, family)
  out <- switch(family,
    Gaussian = pnorm((qnorm(v) - th * qnorm(u2)) / sqrt(1 - th^2)),
    Clayton  = exp(-(th + 1) * log(u2) - (1 + 1 / th) * log(v^(-th) + u2^(-th) - 1)),
    Gumbel   = {
      x <- -log(v); yy <- -log(u2); S <- x^th + yy^th
      exp(-S^(1 / th) + (1 / th - 1) * log(S) + (th - 1) * log(yy) - log(u2))
    },
    Frank    = {
      ## h = e^{-th u2} / ( (e^{-th} - e^{-th u1}) / (e^{-th u1} - 1) + e^{-th u2} ),
      ## a form that avoids the cancellation of the textbook expression for
      ## large |th|
      R <- (exp(-th) - exp(-th * v)) / expm1(-th * v)
      exp(-th * u2) / (R + exp(-th * u2))
    },
    Independence = v
  )
  out <- pmin(pmax(out, 0), 1)           # rounding at extreme theta
  out[u1 <= 0] <- 0; out[u1 >= 1] <- 1
  out[!is.finite(out)] <- NaN
  out
}

## inverse of h in its first argument: u1 such that h(u1 | u2; th) = v
hinv <- function(v, u2, th, family) {
  u2 <- clamp01(u2); v <- clamp01(v)
  th <- clamp_theta(th, family)
  out <- switch(family,
    Gaussian = pnorm(qnorm(v) * sqrt(1 - th^2) + th * qnorm(u2)),
    Clayton  = {
      A <- exp(-th / (1 + th) * (log(v) + (th + 1) * log(u2))) - u2^(-th) + 1
      exp(-log(pmax(A, 1e-300)) / th)
    },
    Frank    = {
      R <- exp(-th * u2) * (1 / v - 1)
      -log((exp(-th) + R) / (R + 1)) / th
    },
    Gumbel   = hinv_bisect(v, u2, th, family),
    Independence = v
  )
  clamp01(out)
}

## vectorised bisection for families without a closed-form inverse; h is
## increasing in u1, and 40 halvings give an absolute error below 1e-12
hinv_bisect <- function(v, u2, th, family, iter = 40L) {
  lo <- rep(0, length(v)); hi <- rep(1, length(v))
  for (k in seq_len(iter)) {
    mid <- (lo + hi) / 2
    up  <- hfunc(mid, u2, th, family) < v
    lo[up] <- mid[up]; hi[!up] <- mid[!up]
  }
  (lo + hi) / 2
}

## ---------------------------------------------------------------- likelihood
## The pieces that the margins contribute at a given W; a block update need
## only recompute those it changes.
marg_pieces <- function(q, y, W, p) {
  lam <- exp(p[["mud"]] + p[["lamd"]] * W)
  my  <- p[["muc"]] + p[["lamc"]] * W
  sig <- sqrt(p[["ssq"]])
  list(u1  = ppois(q, lam), u1m = ppois(q - 1, lam),
       u2  = pnorm(y, my, sig), lf2 = dnorm(y, my, sig, log = TRUE),
       lf1 = dpois(q, lam, log = TRUE))
}

## log p(q_i, y_i | W_i, eta) from cached pieces and copula parameters
loglik_obs <- function(pc, th, family) {
  d <- hfunc(pc$u1, pc$u2, th, family) - hfunc(pc$u1m, pc$u2, th, family)
  out <- log(pmax(d, 1e-300)) + pc$lf2
  out[!is.finite(out)] <- -Inf
  out
}

## the three densities the criteria score, at the sampled W (as in the paper):
##   joint = p(q,y|W), d_g_c = p(q|y,W) = joint / f_c(y|W),
##   c_g_d = p(y|q,W) = joint / f_d(q|W)
logdens_cond_W <- function(q, y, W, p, xi, family, gam = 0.1) {
  pc <- marg_pieces(q, y, W, p)
  th <- theta_of_w(W, xi, family, gam)
  lj <- loglik_obs(pc, th, family)
  cbind(joint = lj, d_g_c = lj - pc$lf2, c_g_d = lj - pc$lf1)
}

## ---------------------------------------------------------------- W integrated out
## p(q_i, y_i | eta) = int p(q_i, y_i | w, eta) phi(w - beta'x_i) dw and the
## corresponding margins, by Simpson quadrature on w = beta'x_i + t, t in
## [-8, 8], step 0.02. Gauss-Hermite rules against the N(beta'x_i, 1) weight
## are not adequate here: for a large count the Poisson factor is sharply
## peaked far from the prior mean, and a 300-node Hermite rule still errs by
## 0.01 in the log density of such an observation, whereas the fixed grid is
## accurate to 1e-6 for every family (test_mixed.R). The margin of y is
## available in closed form,
##   p(y_i | eta) = N(y_i; mu_c + lam_c beta'x_i, sigma^2 + lam_c^2).
simpson_nodes <- function(half = 8, step = 0.02) {
  t <- seq(-half, half, by = step); k <- length(t)
  w <- c(1, rep(c(4, 2), length.out = k - 2), 1) / 3 * step * dnorm(t)
  list(z = t, w = w)
}

## returns an n x 3 matrix of log densities: joint, d_g_c = p(q|y), c_g_d = p(y|q)
logdens_marg_W <- function(q, y, m, p, xi, family, gam = 0.1, gh = simpson_nodes()) {
  n <- length(q)
  acc_j <- numeric(n); acc_q <- numeric(n)
  for (k in seq_along(gh$z)) {
    W  <- m + gh$z[k]
    pc <- marg_pieces(q, y, W, p)
    th <- theta_of_w(W, xi, family, gam)
    acc_j <- acc_j + gh$w[k] * exp(loglik_obs(pc, th, family))
    acc_q <- acc_q + gh$w[k] * exp(pc$lf1)
  }
  lj <- log(pmax(acc_j, 1e-300))
  lq <- log(pmax(acc_q, 1e-300))
  ly <- dnorm(y, p[["muc"]] + p[["lamc"]] * m,
              sqrt(p[["ssq"]] + p[["lamc"]]^2), log = TRUE)
  cbind(joint = lj, d_g_c = lj - ly, c_g_d = lj - lq)
}

## ---------------------------------------------------------------- Welford
## Online mean of p, mean of log p and variance of log p for each observation,
## from which the WAIC-type criteria follow (Section 3.6 of the paper).
crit_state <- function(n, k = 3) list(mp = matrix(0, n, k), ml = matrix(0, n, k),
                                       vl = matrix(0, n, k), cnt = 0)
crit_update <- function(s, L) {                 # L is n x k of log densities
  s$cnt <- s$cnt + 1
  s$mp <- s$mp + (exp(L) - s$mp) / s$cnt
  dl <- L - s$ml
  s$ml <- s$ml + dl / s$cnt
  s$vl <- s$vl + dl * (L - s$ml)
  s
}
## WAIC = -2 lppd + 2 p_waic, per column
crit_waic <- function(s) {
  lppd <- colSums(log(pmax(s$mp, 1e-300)))
  pw   <- colSums(s$vl) / s$cnt
  rbind(lppd = lppd, p_waic = pw, WAIC = -2 * lppd + 2 * pw)
}

## ---------------------------------------------------------------- simulation
## Draw (q, y) given W from the model, using the conditional (Rosenblatt)
## construction: U2 ~ Unif, U1 = hinv(V | U2) with V ~ Unif.
simulate_mixed <- function(W, p, xi, family, gam = 0.1) {
  n <- length(W)
  th <- theta_of_w(W, xi, family, gam)
  u2 <- runif(n); u1 <- hinv(runif(n), u2, th, family)
  q  <- qpois(u1, exp(p[["mud"]] + p[["lamd"]] * W))
  y  <- qnorm(u2, p[["muc"]] + p[["lamc"]] * W, sqrt(p[["ssq"]]))
  list(q = q, y = y, u1 = u1, u2 = u2, theta = th)
}
