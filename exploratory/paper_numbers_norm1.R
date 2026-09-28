## Remaining numbers quoted in the kidney sections of the paper, recomputed on
## the constrained (||beta|| = 1) fits.
source("R/lvc_gamma_gauss.R")

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; n <- length(q)
fams <- c("Gaussian", "Clayton", "Frank", "Gumbel", "Joe", "SurvClayton",
          "Independence")

load_fam <- function(f) lapply(1:2, function(k)
  readRDS(sprintf("output/fits/%s_chain%d.rds", f, k)))

## ---- per-family alpha, lamq, cor(log alpha, xi0) ----------------------------
cat("=== per-family posterior means ===\n")
for (f in fams) {
  fits <- load_fam(f)
  S <- do.call(cbind, lapply(fits, function(fit)
    fit$PAR[, (fit$burn + 1):fit$M]))
  co <- if (f == "Independence") NA else cor(S["lalpha", ], S["xi0", ])
  cat(sprintf("%-13s alpha %6.1f   lamq %.4f   cor(la,xi0) %6.2f\n",
              f, mean(exp(S["lalpha", ])), mean(S["lamq", ]), co))
}

## ---- conditional WAIC penalty (pwaic_c) -------------------------------------
cat("\n=== conditional pWAIC ===\n")
for (f in fams) {
  fits <- load_fam(f)
  cnt <- 0; mlc <- rep(0, n); vlc <- rep(0, n)
  for (fit in fits) {
    nd <- ncol(fit$Wsamp); thin <- (fit$M - fit$burn) / nd
    didx <- fit$burn + seq_len(nd) * thin
    for (j in seq_len(nd)) {
      p <- fit$PAR[1:8, didx[j]]; names(p) <- rownames(fit$PAR)[1:8]
      lc <- lvc_logdens(q, y, fit$Wsamp[, j], p, f, 0.1)
      cnt <- cnt + 1
      d <- lc - mlc; mlc <- mlc + d / cnt; vlc <- vlc + d * (lc - mlc)
    }
  }
  cat(sprintf("%-13s pwaic_c %8.1f\n", f, sum(vlc / cnt)))
}
res <- readRDS("results.rds")
cat("\nmarginal pWAIC range:",
    paste(round(range(res$crit$pwaic_m), 1), collapse = " to "), "\n")

## ---- tau at the 5th/95th percentile of W under Frank ------------------------
fits <- load_fam("Frank")
S <- do.call(cbind, lapply(fits, function(fit)
  fit$PAR[, (fit$burn + 1):fit$M]))
Wall <- do.call(cbind, lapply(fits, function(fit) fit$Wsamp))
w05 <- quantile(Wall, .05); w95 <- quantile(Wall, .95)
debye1 <- function(a) if (abs(a) < 1e-8) 1 else
  integrate(function(t) t / (exp(t) - 1), 0, a, rel.tol = 1e-8)$value / a
tau_frank <- function(th) 1 - 4 / th + 4 / th^2 * (th * sapply(th, debye1))
idx <- round(seq(1, ncol(S), length.out = 2000))
t05 <- sapply(idx, function(j)
  tau_frank(theta_of_w(w05, c(S["xi0", j], S["xi1", j]), "Frank", 0.1)))
t95 <- sapply(idx, function(j)
  tau_frank(theta_of_w(w95, c(S["xi0", j], S["xi1", j]), "Frank", 0.1)))
cat(sprintf("\ntau at 5th pct of W:  mean %.3f\ntau at 95th pct of W: mean %.3f\n",
            mean(t05), mean(t95)))

## ---- W vs age and eGFR ------------------------------------------------------
Wm <- rowMeans(sapply(fits, function(fit) fit$Wmean))
cat(sprintf("\nSpearman(W, age)  = %.3f\n", cor(Wm, NH$dat$age, method = "spearman")))
cat(sprintf("Spearman(W, eGFR) = %.3f\n",
            cor(Wm, NH$dat$egfr, method = "spearman", use = "complete.obs")))
