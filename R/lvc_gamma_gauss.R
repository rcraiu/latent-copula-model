## Gibbs sampler for the gamma-Gaussian latent variable copula model
## (the continuous-continuous canonical model of Section 2.1):
##
##   Q_i | W_i ~ Gamma(alpha, alpha * exp(-(mu_q + lam_q W_i)))   [mean e^{mu_q+lam_q W_i}]
##   Y_i | W_i ~ N(mu_y + lam_y W_i, sigma^2)
##   (Q_i,Y_i) | W_i ~ C(F_q, F_y | theta(W_i, xi))
##   W_i ~ N(beta' x_i, s2W),  ||beta|| = 1
##
## No data augmentation is required: the joint density is available in closed
## form, so D_mis = W only and the Z step of the mixed-response sampler drops out.
##
## Adaptation follows Andrieu & Thoms (2008): componentwise adaptive scaling on
## top of an empirical-covariance proposal for the parameter blocks (their
## Algorithm 6), and global adaptive scaling for each W_i (their Algorithm 4),
## both targeting an acceptance rate of 0.44.

source("R/copula_densities.R")

TARGET_ACC <- 0.44

## ---------------------------------------------------------------- densities
## The three pieces of the joint density are computed separately so that a
## block update need only recompute the pieces it actually changes: the
## q-margin block leaves (lf2, u2) fixed, the y-margin block leaves (lf1, u1)
## fixed, and the copula block leaves both margins fixed.
comp_q <- function(q, W, p) {
  alpha <- exp(p[["lalpha"]])
  ## Parameterised by scale rather than rate so that this agrees bit for bit
  ## with lvc_sampler.cpp. The two are reciprocals, but 1/(alpha/b) and b/alpha
  ## differ in the last bit, and for large shape in the tail of pgamma that
  ## difference is enough to flip a Metropolis decision and send the two
  ## implementations down separate paths.
  scale <- exp(p[["muq"]] + p[["lamq"]] * W) / alpha
  list(lf = dgamma(q, shape = alpha, scale = scale, log = TRUE),
       u  = pgamma(q, shape = alpha, scale = scale))
}
comp_y <- function(y, W, p) {
  my <- p[["muy"]] + p[["lamy"]] * W; sig <- exp(p[["lsigma"]])
  list(lf = dnorm(y, mean = my, sd = sig, log = TRUE),
       u  = pnorm(y, mean = my, sd = sig))
}
comp_th <- function(W, p, family, gscale)
  theta_of_w(W, c(p[["xi0"]], p[["xi1"]]), family, gscale)

## Per-observation log density from cached pieces.
##   "joint" = f(q,y|w),  "q_g_y" = f(q|y,w),  "y_g_q" = f(y|q,w)
logdens_from <- function(cq, cy, th, family, what = "joint") {
  lc <- log_copula_density(cq$u, cy$u, th, family)
  out <- switch(what,
    joint = cq$lf + cy$lf + lc,
    q_g_y = cq$lf + lc,
    y_g_q = cy$lf + lc
  )
  out[!is.finite(out)] <- -Inf
  out
}

## Convenience wrapper (used by the criteria and by external callers).
lvc_logdens <- function(q, y, W, p, family, gscale = 0.1, what = "joint")
  logdens_from(comp_q(q, W, p), comp_y(y, W, p),
               comp_th(W, p, family, gscale), family, what)

## ---------------------------------------------------------------- Welford
new_state <- function(n) list(mp = rep(0, n), ml = rep(0, n), vl = rep(0, n), cnt = 0)
upd_state <- function(s, p_i, logp_i) {
  s$cnt <- s$cnt + 1
  s$mp <- s$mp + (p_i - s$mp) / s$cnt
  dl <- logp_i - s$ml
  s$ml <- s$ml + dl / s$cnt
  s$vl <- s$vl + dl * (logp_i - s$ml)
  s
}
## WAIC-type criterion from a state tracking p and log p
waic_from <- function(s) -2 * sum(log(pmax(s$mp, 1e-300))) + 2 * sum(s$vl / s$cnt)
## CVML-type criterion from a state tracking 1/p and -log p
cvml_from <- function(s) -sum(log(pmax(s$mp, 1e-300)))

## ---------------------------------------------------------------- sampler
run_lvc <- function(q, y, X, family = "Gaussian", M = 20000, burn = NULL,
                    gscale = 0.1, seed = 1729, verbose = 2000,
                    prior = list(s2_mu = 100, s2_lam = 100, s2_xi = 100,
                                 a_alpha = 0.01, b_alpha = 0.01,
                                 a_sig = 0.01, b_sig = 0.01, s2_beta = 1000,
                                 a_w = 0.01, b_w = 0.01),
                    track_criteria = TRUE, thin_W = 0, beta_norm1 = TRUE,
                    wvar_free = TRUE, px_da = TRUE) {
  ## wvar_free frees the conditional variance of W given x, replacing the
  ## normalisation var(W | x) = 1 by ||beta|| = 1; exactly one of the two is
  ## needed, so the sphere constraint is required when the variance is free.
  if (wvar_free && !beta_norm1)
    stop("wvar_free = TRUE requires beta_norm1 = TRUE")
  ## px_da is the parameter-expanded Gibbs step of Xu, Craiu, Sun & Paterson
  ## (2016, JCGS), eq. (3.2). Two reductions are available:
  ##   * identified by var(W|x) = 1 (beta free): draw unrestricted beta* and
  ##     residual scale psi^2, then reduce by psi (the paper's maps).
  ##   * identified by ||beta|| = 1 (s2W free): the same expanded draws, then
  ##     a Metropolis group move on the scale orbit of the identified target.
  if (px_da && beta_norm1 && !wvar_free)
    stop("px_da with beta_norm1 = TRUE requires wvar_free = TRUE")
  if (is.null(prior$a_w)) prior$a_w <- 0.01
  if (is.null(prior$b_w)) prior$b_w <- 0.01
  set.seed(seed)
  n <- length(q); p_dim <- ncol(X)
  if (is.null(burn)) burn <- floor(M / 2)
  has_xi <- family != "Independence"

  pnames <- c("muq", "lamq", "lalpha", "muy", "lamy", "lsigma", "xi0", "xi1")
  PAR  <- matrix(NA_real_, nrow = length(pnames) + p_dim + wvar_free, ncol = M)
  rownames(PAR) <- c(pnames, paste0("beta", 1:p_dim), if (wvar_free) "s2W")

  ## ---- log posterior contributions ------------------------------------
  logprior <- function(p) {
    lp <- dnorm(p["muq"], 0, sqrt(prior$s2_mu), log = TRUE) +
      dnorm(p["lamq"], 0, sqrt(prior$s2_lam), log = TRUE) +
      dnorm(p["muy"], 0, sqrt(prior$s2_mu), log = TRUE) +
      dnorm(p["lamy"], 0, sqrt(prior$s2_lam), log = TRUE) +
      ## alpha ~ Gamma(a,b) with log-scale parametrisation (Jacobian = alpha)
      dgamma(exp(p["lalpha"]), prior$a_alpha, prior$b_alpha, log = TRUE) + p["lalpha"] +
      ## sigma^2 ~ InvGamma(a,b), sampled as log sigma (Jacobian = 2 sigma^2)
      (-prior$a_sig - 1) * (2 * p["lsigma"]) - prior$b_sig / exp(2 * p["lsigma"]) +
      2 * p["lsigma"]
    if (has_xi) lp <- lp + sum(dnorm(c(p["xi0"], p["xi1"]), 0, sqrt(prior$s2_xi), log = TRUE))
    lp
  }
  ## sign constraint for identifiability: lam_y >= 0 fixes the orientation of W
  feasible <- function(p) p[["lamy"]] >= 0

  ll <- function(p, W) sum(lvc_logdens(q, y, W, p, family, gscale))

  ## ---- initialisation (Section 3.1) -----------------------------------
  beta <- as.numeric(solve(crossprod(X), crossprod(X, y)))
  beta <- beta / sqrt(sum(beta^2))          # W has unit conditional variance
  W    <- as.numeric(X %*% beta)

  fit_q <- glm(q ~ W, family = Gamma(link = "log"))
  muq0  <- coef(fit_q)[1]; lamq0 <- coef(fit_q)[2]
  al0   <- 1 / max(summary(fit_q)$dispersion, 1e-3)
  fit_y <- lm(y ~ W)
  muy0  <- coef(fit_y)[1]; lamy0 <- coef(fit_y)[2]
  sig0  <- sqrt(mean(resid(fit_y)^2))

  if (lamy0 < 0) {                          # enforce the sign constraint
    beta <- -beta; W <- -W
    lamq0 <- -lamq0; lamy0 <- -lamy0
  }

  p <- c(muq = as.numeric(muq0), lamq = as.numeric(lamq0), lalpha = log(al0),
         muy = as.numeric(muy0), lamy = as.numeric(lamy0), lsigma = log(sig0),
         xi0 = 0, xi1 = 0)

  if (has_xi) {                             # initialise xi by profile MLE
    obj <- function(v) { pp <- p; pp["xi0"] <- v[1]; pp["xi1"] <- v[2]; -ll(pp, W) }
    op <- try(optim(c(0, 0), obj, method = "Nelder-Mead",
                    control = list(maxit = 400)), silent = TRUE)
    if (!inherits(op, "try-error") && is.finite(op$value)) p[c("xi0", "xi1")] <- op$par
  }

  ## ---- adaptation bookkeeping -----------------------------------------
  blocks <- list(marg_q = c("muq", "lamq", "lalpha"),
                 marg_y = c("muy", "lamy", "lsigma"))
  if (has_xi) blocks$copula <- c("xi0", "xi1")

  ad <- lapply(blocks, function(b) list(
    mu = p[b], Sig = diag(length(b)) * 0.01,
    lam = rep(2.38^2 / length(b), length(b)), acc = rep(0, length(b))
  ))
  lamW <- rep(1, n); muW <- W; varW <- rep(1, n); accW <- rep(0, n)
  gam_t <- function(t) t^(-0.5)

  betapost_prec <- crossprod(X) + diag(1 / prior$s2_beta, p_dim)
  betapost_cov  <- solve(betapost_prec)
  betapost_chol <- chol(betapost_cov)
  ## constrained beta update (||beta|| = 1): projected random-walk Metropolis
  ## on the unit sphere, uniform prior on the sphere, adaptive step size
  sbeta <- 0.1; acc_beta <- 0
  s2W <- 1                                  # conditional variance of W given x
  s_px <- 0.3; acc_px <- 0                  # group-move step on log c

  ## ---- criteria states -------------------------------------------------
  if (track_criteria) {
    st <- list(waic = new_state(n), cwaic_qy = new_state(n), cwaic_yq = new_state(n),
               cvml = new_state(n), ccvml_qy = new_state(n), ccvml_yq = new_state(n))
  }
  ## posterior mean of W, plus thinned draws for post hoc summaries of theta(W)
  Wsum <- rep(0, n); Wcnt <- 0
  Wsamp <- if (thin_W > 0) matrix(NA_real_, n, floor((M - burn) / thin_W)) else NULL
  wj <- 0L

  ## cached density pieces at the current state
  cq <- comp_q(q, W, p); cy <- comp_y(y, W, p); th <- comp_th(W, p, family, gscale)
  cur_ll <- sum(logdens_from(cq, cy, th, family)); cur_lpr <- logprior(p)

  ## which pieces does each block change?
  block_kind <- c(marg_q = "q", marg_y = "y", copula = "th")

  for (t in 1:M) {
    ## ---- Steps 1-3: parameter blocks, componentwise adaptive scaling ---
    for (bi in seq_along(blocks)) {
      b    <- blocks[[bi]]
      kind <- block_kind[[names(blocks)[bi]]]
      sdb  <- sqrt(pmax(diag(ad[[bi]]$Sig), 1e-12))
      for (j in seq_along(b)) {
        prop <- p
        prop[b[j]] <- p[[b[j]]] + sqrt(ad[[bi]]$lam[j]) * sdb[j] * rnorm(1)
        if (!feasible(prop)) {
          a <- 0
        } else {
          nq <- cq; ny <- cy; nth <- th
          if (kind == "q")  nq  <- comp_q(q, W, prop)
          if (kind == "y")  ny  <- comp_y(y, W, prop)
          if (kind == "th") nth <- comp_th(W, prop, family, gscale)
          prop_ll  <- sum(logdens_from(nq, ny, nth, family))
          prop_lpr <- logprior(prop)
          lr <- (prop_ll + prop_lpr) - (cur_ll + cur_lpr)
          a  <- if (is.finite(lr)) min(1, exp(lr)) else 0
          if (runif(1) < a) {
            p <- prop; cur_ll <- prop_ll; cur_lpr <- prop_lpr
            cq <- nq; cy <- ny; th <- nth
            ad[[bi]]$acc[j] <- ad[[bi]]$acc[j] + 1
          }
        }
        ad[[bi]]$lam[j] <- exp(log(ad[[bi]]$lam[j]) + gam_t(t) * (a - TARGET_ACC))
      }
      mu_old <- ad[[bi]]$mu
      ad[[bi]]$mu  <- mu_old + gam_t(t) * (p[b] - mu_old)
      ad[[bi]]$Sig <- ad[[bi]]$Sig + gam_t(t) *
        (tcrossprod(p[b] - mu_old) - ad[[bi]]$Sig)
    }

    ## ---- Step 4: beta, and (if freed) the residual scale of W ----------
    if (px_da && !beta_norm1) {
      ## PX-HC, Xu et al. (2016, eq. 3.2). Working model W ~ N(X beta*, psi^2)
      ## with beta* unrestricted. Reduce by psi:
      ##   beta = beta*/psi,  W = W*/psi,  (lamq, lamy, xi1) = psi * (.*).
      ## After the reduction, var(W | x) = 1 and the scale ridge is a single
      ## Gibbs move in the working parameter.
      prec <- crossprod(X) / s2W + diag(1 / prior$s2_beta, p_dim)
      covb <- solve(prec)
      mb   <- as.numeric(covb %*% (crossprod(X, W) / s2W))
      bstar <- mb + as.numeric(t(chol(covb)) %*% rnorm(p_dim))
      SS   <- sum((W - as.numeric(X %*% bstar))^2)
      psi2 <- 1 / rgamma(1, prior$a_w + n / 2, rate = prior$b_w + SS / 2)
      psi  <- sqrt(psi2)
      beta <- bstar / psi
      W    <- W / psi
      p["lamq"] <- p[["lamq"]] * psi
      p["lamy"] <- p[["lamy"]] * psi
      if (has_xi) p["xi1"] <- p[["xi1"]] * psi
      s2W  <- 1
      muW  <- muW / psi
      varW <- varW / psi^2
      Xb   <- as.numeric(X %*% beta)
      cur_lpr <- logprior(p)
      cq <- comp_q(q, W, p); cy <- comp_y(y, W, p)
      th <- comp_th(W, p, family, gscale)
      cur_ll <- sum(logdens_from(cq, cy, th, family))
    } else if (beta_norm1) {
      ## Metropolis on the unit sphere: propose beta + s z projected back to
      ## the sphere. The proposal density depends only on the angle between
      ## the current and proposed points, hence is symmetric, and the target
      ## on the sphere is exp(-||W - X beta||^2 / (2 s2W)) (flat spherical prior).
      bp <- beta + sbeta * rnorm(p_dim)
      bp <- bp / sqrt(sum(bp^2))
      lr_b <- 0.5 * (sum((W - as.numeric(X %*% beta))^2) -
                     sum((W - as.numeric(X %*% bp))^2)) / s2W
      a_b <- if (is.finite(lr_b)) min(1, exp(lr_b)) else 0
      if (runif(1) < a_b) { beta <- bp; acc_beta <- acc_beta + 1 }
      sbeta <- exp(log(sbeta) + gam_t(t) * (a_b - 0.234))
      Xb <- as.numeric(X %*% beta)
      if (wvar_free) {
        SS  <- sum((W - Xb)^2)
        s2W <- 1 / rgamma(1, prior$a_w + n / 2, rate = prior$b_w + SS / 2)
      }
    } else {
      ## exact conjugate draw, unit residual variance
      bm   <- betapost_cov %*% crossprod(X, W)
      beta <- as.numeric(bm + t(betapost_chol) %*% rnorm(p_dim))
      Xb   <- as.numeric(X %*% beta)
    }

    ## ---- Step 5: W, componentwise adaptive random walk ----------------
    lpW_cur <- logdens_from(cq, cy, th, family) + dnorm(W, Xb, sqrt(s2W), log = TRUE)
    Wp  <- W + sqrt(lamW * pmax(varW, 1e-8)) * rnorm(n)
    pq  <- comp_q(q, Wp, p); py <- comp_y(y, Wp, p)
    pth <- comp_th(Wp, p, family, gscale)
    lpW_prop <- logdens_from(pq, py, pth, family) + dnorm(Wp, Xb, sqrt(s2W), log = TRUE)
    lr  <- lpW_prop - lpW_cur
    aW  <- exp(pmin(lr, 0))
    aW[!is.finite(aW)] <- 0
    take <- runif(n) < aW
    W <- ifelse(take, Wp, W)
    ## splice the accepted components into the cache
    cq$lf <- ifelse(take, pq$lf, cq$lf); cq$u <- ifelse(take, pq$u, cq$u)
    cy$lf <- ifelse(take, py$lf, cy$lf); cy$u <- ifelse(take, py$u, cy$u)
    th    <- if (length(th) == 1L) th else ifelse(take, pth, th)
    accW <- accW + take
    lamW <- exp(log(lamW) + gam_t(t) * (aW - TARGET_ACC))
    muW_old <- muW
    muW  <- muW_old + gam_t(t) * (W - muW_old)
    varW <- varW + gam_t(t) * ((W - muW_old)^2 - varW)
    cur_ll <- sum(logdens_from(cq, cy, th, family))

    ## ---- Step 6: PX group move on the identified (||beta||=1, s2W free) orbit
    ## (W, lamq, lamy, xi1, s2W) -> (c W, lamq/c, lamy/c, xi1/c, c^2 s2W).
    ## The complete-data likelihood is invariant; the Metropolis ratio is the
    ## W-prior, the parameter priors, the IG prior on s2W, and the Jacobian
    ## c^(n - k + 2), with a symmetric proposal on log c.
    if (px_da && beta_norm1) {
      lc <- s_px * rnorm(1); cc <- exp(lc)
      pn <- p
      pn["lamq"] <- p[["lamq"]] / cc
      pn["lamy"] <- p[["lamy"]] / cc
      if (has_xi) pn["xi1"] <- p[["xi1"]] / cc
      s2n <- cc^2 * s2W
      k <- 2 + has_xi
      lr_px <- sum(dnorm(cc * W, Xb, sqrt(s2n), log = TRUE) -
                   dnorm(W,      Xb, sqrt(s2W), log = TRUE)) +
        (logprior(pn) - logprior(p)) +
        (-(prior$a_w + 1) * 2 * lc - prior$b_w * (1 / s2n - 1 / s2W)) +
        (n - k + 2) * lc
      a_px <- if (is.finite(lr_px)) min(1, exp(lr_px)) else 0
      if (runif(1) < a_px) {
        p <- pn; W <- cc * W; s2W <- s2n
        muW <- cc * muW; varW <- cc^2 * varW
        cur_lpr <- logprior(p)
        cq <- comp_q(q, W, p); cy <- comp_y(y, W, p)
        th <- comp_th(W, p, family, gscale)
        cur_ll <- sum(logdens_from(cq, cy, th, family))
        acc_px <- acc_px + 1
      }
      s_px <- exp(log(s_px) + gam_t(t) * (a_px - TARGET_ACC))
    }

    PAR[, t] <- c(p, beta, if (wvar_free) s2W)

    ## ---- criteria, online after burn-in -------------------------------
    if (track_criteria && t > burn) {
      lj <- logdens_from(cq, cy, th, family, "joint")
      lq <- logdens_from(cq, cy, th, family, "q_g_y")
      ly <- logdens_from(cq, cy, th, family, "y_g_q")
      st$waic     <- upd_state(st$waic,     exp(lj), lj)
      st$cwaic_qy <- upd_state(st$cwaic_qy, exp(lq), lq)
      st$cwaic_yq <- upd_state(st$cwaic_yq, exp(ly), ly)
      st$cvml     <- upd_state(st$cvml,     exp(-lj), -lj)
      st$ccvml_qy <- upd_state(st$ccvml_qy, exp(-lq), -lq)
      st$ccvml_yq <- upd_state(st$ccvml_yq, exp(-ly), -ly)
    }
    if (t > burn) {
      Wsum <- Wsum + W; Wcnt <- Wcnt + 1
      if (thin_W > 0 && (t - burn) %% thin_W == 0 && wj < ncol(Wsamp)) {
        wj <- wj + 1L; Wsamp[, wj] <- W
      }
    }

    if (verbose && t %% verbose == 0)
      cat(sprintf("  [%s] %6d/%d  loglik %.1f\n", family, t, M, cur_ll))
  }

  out <- list(PAR = PAR, family = family, M = M, burn = burn, W = W,
              Wmean = Wsum / max(Wcnt, 1), Wsamp = Wsamp,
              acc = list(blocks = lapply(ad, function(a) a$acc / M),
                         W = mean(accW) / M,
                         beta = if (beta_norm1) acc_beta / M else NA,
                         px = if (px_da && beta_norm1) acc_px / M else NA))
  if (track_criteria) {
    out$criteria <- c(
      WAIC      = waic_from(st$waic),
      CWAIC_qy  = waic_from(st$cwaic_qy),
      CWAIC_yq  = waic_from(st$cwaic_yq),
      CVML      = cvml_from(st$cvml),
      CCVML_qy  = cvml_from(st$ccvml_qy),
      CCVML_yq  = cvml_from(st$ccvml_yq)
    )
  }
  out
}

## posterior summary helper
lvc_summary <- function(fit, probs = c(0.025, 0.5, 0.975)) {
  S <- fit$PAR[, (fit$burn + 1):fit$M, drop = FALSE]
  if (fit$family == "Independence") S <- S[!rownames(S) %in% c("xi0", "xi1"), ]
  tab <- t(apply(S, 1, function(v) c(mean = mean(v), sd = sd(v), quantile(v, probs))))
  ## report alpha and sigma^2 on their natural scales as well
  extra <- rbind(
    alpha = c(mean(exp(S["lalpha", ])), sd(exp(S["lalpha", ])),
              quantile(exp(S["lalpha", ]), probs)),
    sigma2 = c(mean(exp(2 * S["lsigma", ])), sd(exp(2 * S["lsigma", ])),
               quantile(exp(2 * S["lsigma", ]), probs))
  )
  colnames(extra) <- colnames(tab)
  rbind(tab, extra)
}
