## Gibbs sampler for the mixed-response latent variable copula model of
## Section 3 of the paper, for the Gaussian, Clayton, Frank, Gumbel and
## independence copulas. Run from the repository root.
##
## One sweep updates, in this order,
##   1. (mu_d, lam_d)            RW Metropolis on the observed-data likelihood
##   2. (mu_c, lam_c, log s2)    RW Metropolis on the observed-data likelihood
##   3. U_i                      the augmented continuous variable behind the
##                               count: U_i | rest is the conditional copula
##                               restricted to (F_d(q_i - 1), F_d(q_i)], drawn
##                               by inverting h(. | F_c(y_i))
##   4. xi                       RW Metropolis on the augmented likelihood
##                               prod_i c(U_i, F_c(y_i) | theta(W_i, xi))
##   5. beta                     conjugate normal draw given W
##   6. W_i                      componentwise RW Metropolis
##   7. scale group move         (W, beta, lam_d, lam_c, xi_1) ->
##                               (cW, c beta, lam_d/c, lam_c/c, xi_1/c), a
##                               parameter-expanded move along the ridge that
##                               trades the scale of W against the loadings; the
##                               observed-data likelihood is invariant and the
##                               ratio involves only the prior of W given beta,
##                               the parameter priors and the Jacobian
## Steps 1, 2 and 6 use the likelihood with U integrated out, and step 3
## redraws U from its full conditional before step 4 uses it; this is a
## partially collapsed Gibbs sampler and leaves the augmented posterior
## invariant. Section 3.4 of the paper describes the same moves with the
## augmented variable on the normal scale, Z_i = Phi^{-1}(U_i).
##
## Adaptation follows Andrieu & Thoms (2008): componentwise adaptive scaling
## with an adapted covariance for the parameter blocks (their Algorithm 6) and
## global adaptive scaling for each W_i (their Algorithm 4), both targeting an
## acceptance rate of 0.44.
##
## The criteria as defined in the paper, which score p(q_i, y_i | eta, W_i) at
## the sampled W_i, are accumulated online during the run. The W-marginalised
## criteria that the paper now reports are computed afterwards from the saved
## parameter draws by criteria.R, since they need a quadrature per draw and
## per observation.

source("mixed-response/mixed_model.R")

TARGET_ACC <- 0.44

run_mixed <- function(q, y, X, family = "Gaussian", M = 100000, burn = NULL,
                      seed = 1729, gam = 0.1,
                      prior = list(s2_beta = 1000, s2_mu = 100, s2_lam = 100,
                                   a_sig = 1e-4, b_sig = 1e-4,
                                   mu_xi = c(0, 0), s2_xi = 100),
                      jitter = 0, track_criteria = TRUE, verbose = 10000) {
  stopifnot(family %in% c("Gaussian", "Clayton", "Frank", "Gumbel", "Independence"))
  set.seed(seed)
  n <- length(q); pd <- ncol(X)
  if (is.null(burn)) burn <- floor(M / 2)
  has_xi <- family != "Independence"

  pnames <- c("mud", "lamd", "muc", "lamc", "ssq", "xi0", "xi1")
  PAR <- matrix(NA_real_, length(pnames) + pd, M)
  rownames(PAR) <- c(paste0("beta", seq_len(pd)), pnames)

  ## ---- priors ------------------------------------------------------------
  logprior <- function(p) {
    lp <- dnorm(p[["mud"]], 0, sqrt(prior$s2_mu), log = TRUE) +
      dnorm(p[["lamd"]], 0, sqrt(prior$s2_lam), log = TRUE) +
      dnorm(p[["muc"]], 0, sqrt(prior$s2_mu), log = TRUE) +
      dnorm(p[["lamc"]], 0, sqrt(prior$s2_lam), log = TRUE) +
      ## sigma^2 ~ InvGamma(a, b), proposed on the log scale: the Jacobian
      ## sigma^2 is included here so that the ratio is taken in log sigma^2
      (-prior$a_sig - 1) * log(p[["ssq"]]) - prior$b_sig / p[["ssq"]] + log(p[["ssq"]])
    if (has_xi)
      lp <- lp + sum(dnorm(c(p[["xi0"]], p[["xi1"]]), prior$mu_xi, sqrt(prior$s2_xi), log = TRUE))
    lp
  }
  feasible <- function(p) p[["lamc"]] >= 0 && p[["ssq"]] > 0   # sign identification

  ## ---- initialisation (Section 3.1) --------------------------------------
  beta <- as.numeric(solve(crossprod(X), crossprod(X, y)))
  W    <- as.numeric(X %*% beta)
  if (jitter > 0) {                        # overdispersed start for later chains
    beta <- beta * exp(jitter * rnorm(1)) + jitter * rnorm(pd)
    W    <- as.numeric(X %*% beta)
  }
  gq <- glm(q ~ W, family = poisson())
  ly <- lm(y ~ W)
  p <- c(mud = unname(coef(gq)[1]), lamd = unname(coef(gq)[2]),
         muc = unname(coef(ly)[1]), lamc = unname(coef(ly)[2]),
         ssq = mean(resid(ly)^2), xi0 = 0, xi1 = 0)
  if (p[["lamc"]] < 0) {                   # orient W so that lam_c >= 0
    beta <- -beta; W <- -W; p["lamd"] <- -p[["lamd"]]; p["lamc"] <- -p[["lamc"]]
  }
  if (jitter > 0) {
    p[c("mud", "lamd", "muc")] <- p[c("mud", "lamd", "muc")] + jitter * rnorm(3)
    p["lamc"] <- abs(p[["lamc"]] + jitter * rnorm(1))
    p["ssq"]  <- p[["ssq"]] * exp(jitter * rnorm(1))
  }
  pc <- marg_pieces(q, y, W, p)
  if (has_xi) {                            # profile optimum for xi given the rest
    obj <- function(v) -sum(loglik_obs(pc, theta_of_w(W, v, family, gam), family))
    op <- try(optim(c(0, 0), obj, method = "Nelder-Mead", control = list(maxit = 500)),
              silent = TRUE)
    if (!inherits(op, "try-error") && is.finite(op$value)) p[c("xi0", "xi1")] <- op$par
    if (jitter > 0) p[c("xi0", "xi1")] <- p[c("xi0", "xi1")] + jitter * rnorm(2)
  }
  xi <- c(p[["xi0"]], p[["xi1"]])
  th <- theta_of_w(W, xi, family, gam)
  U  <- (pc$u1 + pc$u1m) / 2

  ## ---- adaptation bookkeeping ------------------------------------------
  blocks <- list(marg_d = c("mud", "lamd"), marg_c = c("muc", "lamc", "ssq"))
  if (has_xi) blocks$copula <- c("xi0", "xi1")
  ## the variance block is adapted on the log scale
  tolog <- function(b, v) { if ("ssq" %in% b) v["ssq"] <- log(v[["ssq"]]); v }
  ad <- lapply(blocks, function(b) list(
    mu = tolog(b, p[b]), Sig = diag(length(b)) * 0.01,
    lam = rep(1, length(b)), acc = rep(0, length(b))))
  lamW <- rep(0.5, n); muW <- W; varW <- rep(1, n); accW <- rep(0, n)
  lsc <- log(0.05); accS <- 0                # scale group move, adapted on log c
  gam_t <- function(t) t^(-0.6)

  bpost_cov  <- solve(crossprod(X) + diag(1 / prior$s2_beta, pd))
  bpost_chol <- t(chol(bpost_cov))

  cur_ll  <- sum(loglik_obs(pc, th, family))
  cur_lpr <- logprior(p)

  if (track_criteria) st <- crit_state(n)
  Wsum <- rep(0, n)

  for (t in seq_len(M)) {
    ## ---- Steps 1-2: marginal blocks, componentwise adaptive scaling ------
    for (bn in c("marg_d", "marg_c")) {
      b <- blocks[[bn]]
      sdb <- sqrt(pmax(diag(ad[[bn]]$Sig), 1e-10))
      for (j in seq_along(b)) {
        prop <- p
        step <- sqrt(ad[[bn]]$lam[j]) * sdb[j] * rnorm(1)
        if (b[j] == "ssq") prop["ssq"] <- p[["ssq"]] * exp(step) else prop[b[j]] <- p[[b[j]]] + step
        a <- 0
        if (feasible(prop)) {
          npc <- marg_pieces(q, y, W, prop)
          prop_ll <- sum(loglik_obs(npc, th, family)); prop_lpr <- logprior(prop)
          lr <- (prop_ll + prop_lpr) - (cur_ll + cur_lpr)
          a <- if (is.finite(lr)) min(1, exp(lr)) else 0
          if (runif(1) < a) {
            p <- prop; pc <- npc; cur_ll <- prop_ll; cur_lpr <- prop_lpr
            ad[[bn]]$acc[j] <- ad[[bn]]$acc[j] + 1
          }
        }
        ad[[bn]]$lam[j] <- exp(log(ad[[bn]]$lam[j]) + gam_t(t) * (a - TARGET_ACC))
      }
      v <- tolog(b, p[b]); mu_old <- ad[[bn]]$mu
      ad[[bn]]$mu  <- mu_old + gam_t(t) * (v - mu_old)
      ad[[bn]]$Sig <- ad[[bn]]$Sig + gam_t(t) * (tcrossprod(v - mu_old) - ad[[bn]]$Sig)
    }

    if (has_xi) {
      ## ---- Step 3: augmented variable U_i | rest --------------------------
      hlo <- hfunc(pc$u1m, pc$u2, th, family); hhi <- hfunc(pc$u1, pc$u2, th, family)
      V <- hlo + runif(n) * pmax(hhi - hlo, 0)
      U <- hinv(V, pc$u2, th, family)
      U <- pmin(pmax(U, pc$u1m), pc$u1)             # guard against rounding
      U <- clamp01(U)

      ## ---- Step 4: xi | U, W, margins ------------------------------------
      b <- blocks$copula
      cur_lc <- sum(log_copula_density(U, pc$u2, th, family))
      sdb <- sqrt(pmax(diag(ad$copula$Sig), 1e-10))
      for (j in seq_along(b)) {
        prop <- p
        prop[b[j]] <- p[[b[j]]] + sqrt(ad$copula$lam[j]) * sdb[j] * rnorm(1)
        nth <- theta_of_w(W, c(prop[["xi0"]], prop[["xi1"]]), family, gam)
        prop_lc <- sum(log_copula_density(U, pc$u2, nth, family)); prop_lpr <- logprior(prop)
        lr <- (prop_lc + prop_lpr) - (cur_lc + cur_lpr)
        a <- if (is.finite(lr)) min(1, exp(lr)) else 0
        if (runif(1) < a) {
          p <- prop; th <- nth; cur_lc <- prop_lc; cur_lpr <- prop_lpr
          ad$copula$acc[j] <- ad$copula$acc[j] + 1
        }
        ad$copula$lam[j] <- exp(log(ad$copula$lam[j]) + gam_t(t) * (a - TARGET_ACC))
      }
      v <- p[b]; mu_old <- ad$copula$mu
      ad$copula$mu  <- mu_old + gam_t(t) * (v - mu_old)
      ad$copula$Sig <- ad$copula$Sig + gam_t(t) * (tcrossprod(v - mu_old) - ad$copula$Sig)
      xi <- c(p[["xi0"]], p[["xi1"]])
      cur_ll <- sum(loglik_obs(pc, th, family))
    }

    ## ---- Step 5: beta | W, conjugate ---------------------------------------
    beta <- as.numeric(bpost_cov %*% crossprod(X, W) + bpost_chol %*% rnorm(pd))
    Xb   <- as.numeric(X %*% beta)

    ## ---- Step 6: W_i, componentwise adaptive random walk ------------------
    lpW_cur <- loglik_obs(pc, th, family) + dnorm(W, Xb, 1, log = TRUE)
    Wp  <- W + sqrt(lamW * pmax(varW, 1e-8)) * rnorm(n)
    ppc <- marg_pieces(q, y, Wp, p)
    pth <- theta_of_w(Wp, xi, family, gam)
    lpW_prop <- loglik_obs(ppc, pth, family) + dnorm(Wp, Xb, 1, log = TRUE)
    aW <- exp(pmin(lpW_prop - lpW_cur, 0)); aW[!is.finite(aW)] <- 0
    take <- runif(n) < aW
    W <- ifelse(take, Wp, W)
    for (nm in names(pc)) pc[[nm]] <- ifelse(take, ppc[[nm]], pc[[nm]])
    th <- ifelse(take, pth, th)
    accW <- accW + take
    lamW <- exp(log(lamW) + gam_t(t) * (aW - TARGET_ACC))
    muW_old <- muW
    muW  <- muW_old + gam_t(t) * (W - muW_old)
    varW <- varW + gam_t(t) * ((W - muW_old)^2 - varW)
    cur_ll <- sum(loglik_obs(pc, th, family))

    ## ---- Step 7: scale group move -----------------------------------------
    ## c = exp(eps), eps ~ N(0, s^2) symmetric; the complete-data likelihood
    ## is invariant, so the log ratio is the change in log N(W; X beta, 1),
    ## in the priors of beta, lam_d, lam_c and xi_1, plus the Jacobian
    ## (n + pd - k) log c with k the number of loadings that scale.
    eps <- exp(lsc) * rnorm(1); cc <- exp(eps)
    r2  <- sum((W - Xb)^2); k <- if (has_xi) 3 else 2
    lr  <- -(cc^2 - 1) / 2 * r2 - (cc^2 - 1) * sum(beta^2) / (2 * prior$s2_beta) -
      (1 / cc^2 - 1) * (p[["lamd"]]^2 + p[["lamc"]]^2) / (2 * prior$s2_lam) +
      (n + pd - k) * eps
    if (has_xi)
      lr <- lr - ((p[["xi1"]] / cc - prior$mu_xi[2])^2 - (p[["xi1"]] - prior$mu_xi[2])^2) /
        (2 * prior$s2_xi)
    aS <- if (is.finite(lr)) min(1, exp(lr)) else 0
    if (runif(1) < aS) {
      W <- cc * W; beta <- cc * beta; Xb <- cc * Xb
      p["lamd"] <- p[["lamd"]] / cc; p["lamc"] <- p[["lamc"]] / cc
      if (has_xi) { p["xi1"] <- p[["xi1"]] / cc; xi <- c(p[["xi0"]], p[["xi1"]]) }
      muW <- cc * muW; varW <- cc^2 * varW
      pc <- marg_pieces(q, y, W, p); th <- theta_of_w(W, xi, family, gam)
      cur_ll <- sum(loglik_obs(pc, th, family)); cur_lpr <- logprior(p)
      accS <- accS + 1
    }
    lsc <- lsc + gam_t(t) * (aS - TARGET_ACC)

    PAR[, t] <- c(beta, p)

    ## ---- criteria as defined in the paper, at the sampled W ---------------
    if (t > burn) {
      Wsum <- Wsum + W
      if (track_criteria) {
        lj <- loglik_obs(pc, th, family)
        st <- crit_update(st, cbind(lj, lj - pc$lf2, lj - pc$lf1))
      }
    }
    if (verbose && t %% verbose == 0)
      cat(sprintf("  [%s] %7d/%d  loglik %.1f  acc W %.2f  acc scale %.2f\n", family, t, M,
                  cur_ll, mean(accW) / t, accS / t))
  }

  out <- list(PAR = PAR, family = family, M = M, burn = burn, seed = seed, gam = gam,
              prior = prior, Wmean = Wsum / (M - burn),
              acc = list(blocks = lapply(ad, function(a) a$acc / M), W = mean(accW) / M,
                         scale = accS / M))
  if (track_criteria) {
    cw <- crit_waic(st)
    colnames(cw) <- c("joint", "d_g_c", "c_g_d")
    out$criteria_condW <- cw
  }
  out
}

## posterior means and the criteria pooled over several chains
pool_chains <- function(fits) {
  S <- do.call(cbind, lapply(fits, function(f) f$PAR[, (f$burn + 1):f$M, drop = FALSE]))
  list(means = rowMeans(S), draws = S)
}
