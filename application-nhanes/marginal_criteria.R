## Marginal (W integrated out) predictive criteria.
##
## The conditional criteria of the paper evaluate p(q_i,y_i | eta, W_i). Because
## the model carries one latent variable per observation, W_i acts as a free
## per-observation parameter and the conditional density can be made large by
## shrinking the marginal dispersions and letting W_i chase each data point.
## The quantity a new observation is actually drawn from is
##      p(q_i,y_i | eta) = int f(q_i,y_i | w) phi(w - beta'x_i) dw,
## which we evaluate by adaptive Simpson quadrature. The integrand is sharply
## peaked (its width is of order 1/(lam_q sqrt(alpha))), so we first locate its
## mode on a coarse grid and then integrate on a fine grid around the mode.

source("R/lvc_gamma_gauss.R")

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X; n <- length(q)
fams <- c("Gaussian", "Clayton", "Frank", "Gumbel", "Independence")

simpson_w <- function(k) {                       # weights for k odd points
  w <- c(1, rep(c(4, 2), length.out = k - 2), 1)
  w / 3
}

## log integrand, vectorised over the n observations at a single w-offset
log_integrand <- function(w, p, m, family)
  lvc_logdens(q, y, w, p, family, 0.1, "joint") - 0.5 * (w - m)^2 - 0.9189385332

log_marg <- function(p, beta, family, half = 1.6, kf = 81, kc = 49) {
  m <- as.numeric(X %*% beta)
  ## coarse pass: locate the mode of the integrand for each observation
  gc <- seq(-6, 6, length.out = kc)
  best <- matrix(-Inf, n, 1); what <- numeric(n)
  for (t in gc) {
    v <- log_integrand(m + t, p, m, family)
    hit <- v > best[, 1]
    best[hit, 1] <- v[hit]; what[hit] <- t
  }
  ## fine pass: Simpson on [mode - half, mode + half]
  tt <- seq(-half, half, length.out = kf)
  h  <- tt[2] - tt[1]
  sw <- simpson_w(kf)
  mx <- best[, 1]                                # for log-sum-exp stability
  acc <- numeric(n)
  for (j in seq_len(kf)) {
    v <- log_integrand(m + what + tt[j], p, m, family)
    acc <- acc + sw[j] * exp(v - mx)
  }
  mx + log(h * acc)
}

## ---- accuracy check on a subsample, coarse vs fine ------------------------
fit <- readRDS("output/fits/Independence_chain1.rds")
p0 <- fit$PAR[1:8, fit$M]; names(p0) <- rownames(fit$PAR)[1:8]
b0 <- fit$PAR[9:12, fit$M]
a1 <- log_marg(p0, b0, "Independence", half = 1.6, kf = 81)
a2 <- log_marg(p0, b0, "Independence", half = 2.4, kf = 401)
cat("quadrature check (Independence, tightest integrand):\n")
cat("  max |log p difference| between 81 and 401 nodes:",
    signif(max(abs(a1 - a2)), 3), "\n")
cat("  sum difference:", signif(sum(a1) - sum(a2), 4), "on a total of",
    round(sum(a2), 1), "\n\n")

## ---- criteria -------------------------------------------------------------
## Wsamp column j holds W at iteration burn + j*thin_W, so the parameter draw
## that goes with it is PAR[, burn + j*thin_W]. Pairing W draws with the wrong
## parameter draws silently destroys both criteria, so the mapping is built
## from the stored dimensions rather than assumed.
res <- NULL
for (f in fams) {
  fit <- readRDS(sprintf("output/fits/%s_chain%d.rds", f, 1))
  Sp  <- fit$PAR
  nd  <- ncol(fit$Wsamp)
  thin <- (fit$M - fit$burn) / nd
  didx <- fit$burn + seq_len(nd) * thin
  stopifnot(max(didx) <= fit$M, min(didx) > fit$burn)

  cnt <- 0
  mp <- rep(0, n);  ml <- rep(0, n);  vl <- rep(0, n);  mip <- rep(0, n)
  mpc <- rep(0, n); mlc <- rep(0, n); vlc <- rep(0, n); mipc <- rep(0, n)
  for (j in seq_len(nd)) {
    d <- didx[j]
    p <- Sp[1:8, d]; names(p) <- rownames(Sp)[1:8]
    lj <- log_marg(p, Sp[9:12, d], f)                       # W integrated out
    lc <- lvc_logdens(q, y, fit$Wsamp[, j], p, f, 0.1)      # W conditioned on
    cnt <- cnt + 1
    mp  <- mp  + (exp(lj) - mp) / cnt
    mip <- mip + (exp(-lj) - mip) / cnt
    d1 <- lj - ml; ml <- ml + d1 / cnt; vl <- vl + d1 * (lj - ml)
    mpc  <- mpc  + (exp(lc) - mpc) / cnt
    mipc <- mipc + (exp(-lc) - mipc) / cnt
    d2 <- lc - mlc; mlc <- mlc + d2 / cnt; vlc <- vlc + d2 * (lc - mlc)
  }
  lppd <- sum(log(pmax(mp, 1e-300))); pen <- sum(vl / cnt)
  lppdc <- sum(log(pmax(mpc, 1e-300))); penc <- sum(vlc / cnt)
  res <- rbind(res, data.frame(
    family = f, draws = cnt,
    lppd_m = lppd, pwaic_m = pen, WAIC_m = -2 * lppd + 2 * pen,
    CVML_m = -sum(log(mip)),
    lppd_c = lppdc, pwaic_c = penc, WAIC_c = -2 * lppdc + 2 * penc,
    CVML_c = -sum(log(mipc))))
  cat(sprintf("%-13s marg: lppd %8.1f pwaic %6.1f WAIC %9.1f | cond: lppd %8.1f pwaic %7.1f WAIC %9.1f\n",
              f, lppd, pen, -2 * lppd + 2 * pen, lppdc, penc, -2 * lppdc + 2 * penc))
}

cat("\n=== marginal criteria, ranked (lower WAIC, higher CVML is better) ===\n")
rm_ <- res[order(res$WAIC_m), ]
print(data.frame(family = rm_$family, WAIC = round(rm_$WAIC_m, 1),
                 dWAIC = round(rm_$WAIC_m - min(rm_$WAIC_m), 1),
                 CVML = round(rm_$CVML_m, 1), row.names = NULL))

cat("\n=== conditional criteria (as defined in the paper), ranked ===\n")
rc_ <- res[order(res$WAIC_c), ]
print(data.frame(family = rc_$family, WAIC = round(rc_$WAIC_c, 1),
                 dWAIC = round(rc_$WAIC_c - min(rc_$WAIC_c), 1),
                 CVML = round(rc_$CVML_c, 1), row.names = NULL))
saveRDS(res, "output/marginal_criteria.rds")
