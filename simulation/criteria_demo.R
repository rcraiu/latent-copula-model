## Demonstration that the conditional criteria degenerate when both margins
## carry a free dispersion parameter.
##
## Data are simulated from the canonical gamma-Gaussian model with a Gaussian
## copula whose parameter genuinely depends on W, so the correct answer is known:
## the Gaussian copula model should beat the independence model. We fit both and
## score them by (a) the conditional criteria of the paper, which evaluate
## p(q_i,y_i | eta, W_i), and (b) their marginalised counterparts, which
## integrate W out. Only the latter should recover the truth.

source("R/lvc_gamma_gauss.R")

set.seed(20260907)
n <- 1500
p <- 2
X <- cbind(rnorm(n), rnorm(n))
beta_true <- c(0.6, -0.4)
W <- as.numeric(X %*% beta_true) + rnorm(n)

muq_t <- 0.0; lamq_t <- 0.5; alpha_t <- 12
muy_t <- 0.0; lamy_t <- 0.8; sigma_t <- 0.7
xi_t   <- c(6, 3)                    # theta = tanh(0.1*(6 + 3W)), well away from 0

th <- theta_of_w(W, xi_t, "Gaussian", 0.1)
z1 <- rnorm(n); z2 <- th * z1 + sqrt(1 - th^2) * rnorm(n)
u1 <- pnorm(z1); u2 <- pnorm(z2)
q <- qgamma(u1, shape = alpha_t, rate = alpha_t / exp(muq_t + lamq_t * W))
y <- qnorm(u2, mean = muy_t + lamy_t * W, sd = sigma_t)

cat("simulated n =", n, " Kendall tau =", round(cor(q, y, method = "kendall"), 3),
    " theta range", paste(round(range(th), 3), collapse = " to "), "\n\n")

M <- as.integer(Sys.getenv("DEMO_M", "60000"))
fits <- list()
for (f in c("Gaussian", "Independence")) {
  cat("fitting", f, "...\n")
  fits[[f]] <- run_lvc(q, y, X, family = f, M = M, burn = floor(M / 2),
                       gscale = 0.1, seed = 11, verbose = 0,
                       thin_W = floor(M / 2 / 300))
}

## ---- criteria -------------------------------------------------------------
simpson_w <- function(k) c(1, rep(c(4, 2), length.out = k - 2), 1) / 3
log_marg <- function(pr, beta, family, half = 1.8, kf = 81, kc = 49) {
  m <- as.numeric(X %*% beta)
  li <- function(t) lvc_logdens(q, y, m + t, pr, family, 0.1) - 0.5 * t^2 - 0.9189385332
  best <- rep(-Inf, n); what <- numeric(n)
  for (t in seq(-6, 6, length.out = kc)) {
    v <- li(t); h <- v > best; best[h] <- v[h]; what[h] <- t
  }
  tt <- seq(-half, half, length.out = kf); hh <- tt[2] - tt[1]; sw <- simpson_w(kf)
  acc <- numeric(n)
  for (j in seq_len(kf)) acc <- acc + sw[j] * exp(li(what + tt[j]) - best)
  best + log(hh * acc)
}

out <- NULL
for (f in names(fits)) {
  fit <- fits[[f]]; nd <- ncol(fit$Wsamp)
  thin <- (fit$M - fit$burn) / nd; didx <- fit$burn + seq_len(nd) * thin
  cnt <- 0
  mp <- rep(0, n); mlm <- rep(0, n); vlm <- rep(0, n)
  mpc <- rep(0, n); mlc <- rep(0, n); vlc <- rep(0, n)
  for (j in seq_len(nd)) {
    pr <- fit$PAR[1:8, didx[j]]; names(pr) <- rownames(fit$PAR)[1:8]
    lj <- log_marg(pr, fit$PAR[9:(8 + p), didx[j]], f)
    lc <- lvc_logdens(q, y, fit$Wsamp[, j], pr, f, 0.1)
    cnt <- cnt + 1
    mp <- mp + (exp(lj) - mp) / cnt
    d1 <- lj - mlm; mlm <- mlm + d1 / cnt; vlm <- vlm + d1 * (lj - mlm)
    mpc <- mpc + (exp(lc) - mpc) / cnt
    d2 <- lc - mlc; mlc <- mlc + d2 / cnt; vlc <- vlc + d2 * (lc - mlc)
  }
  lp <- sum(log(pmax(mp, 1e-300))); pw <- sum(vlm / cnt)
  lpc <- sum(log(pmax(mpc, 1e-300))); pwc <- sum(vlc / cnt)
  out <- rbind(out, data.frame(family = f,
    alpha = mean(exp(fit$PAR["lalpha", (fit$burn + 1):fit$M])),
    sigma = mean(exp(fit$PAR["lsigma", (fit$burn + 1):fit$M])),
    WAIC_cond = -2 * lpc + 2 * pwc, pwaic_cond = pwc,
    WAIC_marg = -2 * lp + 2 * pw, pwaic_marg = pw))
}

cat("\ntrue alpha =", alpha_t, " true sigma =", sigma_t, "\n\n")
print(out, digits = 5, row.names = FALSE)
cat("\nconditional WAIC selects:", out$family[which.min(out$WAIC_cond)], "\n")
cat("marginal    WAIC selects:", out$family[which.min(out$WAIC_marg)],
    "  (truth: Gaussian)\n")
saveRDS(out, "output/criteria_demo.rds")
