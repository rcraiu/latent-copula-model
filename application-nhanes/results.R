## Final results for the NHANES continuous-continuous application.
##   * convergence diagnostics across the two chains per family
##   * model comparison on both the conditional criteria of the paper and their
##     marginalised counterparts (W integrated out by quadrature)
##   * posterior summaries, and the implied dependence as a function of W
##   * LaTeX tables and figures
##
## Usage: Rscript results.R [fitdir]

source("R/lvc_gamma_gauss.R")
suppressMessages({library(coda); library(ggplot2)})

FITDIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "output/fits"
dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X; n <- length(q)
fams <- c("Gaussian", "Clayton", "Frank", "Gumbel", "Joe", "SurvClayton",
          "Independence")
CHAINS <- 1:2

## Restrict to families with a complete set of chains, so that the script is
## usable while some fits are still running.
have <- vapply(fams, function(f)
  all(file.exists(file.path(FITDIR, sprintf("%s_chain%d.rds", f, CHAINS)))),
  logical(1))
if (any(!have))
  cat("skipping (no fits):", paste(fams[!have], collapse = ", "), "\n")
fams <- fams[have]
stopifnot(length(fams) > 0)

## ---------------------------------------------------------------- criteria --
simpson_w <- function(k) c(1, rep(c(4, 2), length.out = k - 2), 1) / 3

log_marg <- function(p, beta, family, half = 1.6, kf = 81, kc = 49) {
  m <- as.numeric(X %*% beta)
  li <- function(t) lvc_logdens(q, y, m + t, p, family, 0.1) -
                    0.5 * t^2 - 0.9189385332
  best <- rep(-Inf, n); what <- numeric(n)
  for (t in seq(-6, 6, length.out = kc)) {
    v <- li(t); hit <- v > best; best[hit] <- v[hit]; what[hit] <- t
  }
  tt <- seq(-half, half, length.out = kf); h <- tt[2] - tt[1]; sw <- simpson_w(kf)
  acc <- numeric(n)
  for (j in seq_len(kf)) acc <- acc + sw[j] * exp(li(what + tt[j]) - best)
  best + log(h * acc)
}

## ------------------------------------------------------- Kendall's tau map --
debye1 <- function(th) sapply(th, function(a)
  if (abs(a) < 1e-8) 1 else integrate(function(t) t / (exp(t) - 1), 0, a,
                                      rel.tol = 1e-8)$value / a)
tau_of_theta <- function(th, family) switch(family,
  Gaussian     = 2 / pi * asin(pmax(pmin(th, 1), -1)),
  Clayton      = th / (th + 2),
  Gumbel       = 1 - 1 / th,
  Frank        = 1 - 4 / th + 4 / th^2 * (th * debye1(th)),
  Independence = rep(0, length(th)))

## ------------------------------------------------------------ load & pool ---
load_fam <- function(f) {
  fits <- lapply(CHAINS, function(k)
    readRDS(sprintf("%s/%s_chain%d.rds", FITDIR, f, k)))
  post <- lapply(fits, function(fit) fit$PAR[, (fit$burn + 1):fit$M, drop = FALSE])
  list(fits = fits, post = post)
}

pretty <- c(muq = "$\\mu_\\q$", lamq = "$\\lambda_\\q$", alpha = "$\\alpha$",
            muy = "$\\mu_\\y$", lamy = "$\\lambda_\\y$", sigma = "$\\sigma$",
            xi0 = "$\\xi_0$", xi1 = "$\\xi_1$",
            beta1 = "$\\beta_\\text{age}$", beta2 = "$\\beta_\\text{female}$",
            beta3 = "$\\beta_\\text{BMI}$", beta4 = "$\\beta_\\text{SBP}$")

crit <- NULL; summ <- list(); diagn <- NULL
for (f in fams) {
  L <- load_fam(f); fits <- L$fits

  ## --- convergence ---------------------------------------------------------
  ml <- mcmc.list(lapply(L$post, function(S) {
    S2 <- S; S2["lalpha", ] <- exp(S["lalpha", ]); S2["lsigma", ] <- exp(S["lsigma", ])
    rownames(S2)[rownames(S2) %in% c("lalpha", "lsigma")] <- c("alpha", "sigma")
    if (f == "Independence") S2 <- S2[!rownames(S2) %in% c("xi0", "xi1"), ]
    mcmc(t(S2[, seq(1, ncol(S2), by = 20)]))
  }))
  g <- gelman.diag(ml, multivariate = TRUE, autoburnin = FALSE)
  diagn <- rbind(diagn, data.frame(family = f, maxPSRF = max(g$psrf[, 1]),
                                   MPSRF = g$mpsrf, minESS = min(effectiveSize(ml))))

  ## --- criteria, pooling the two chains ------------------------------------
  cnt <- 0
  mp <- rep(0, n); ml_ <- rep(0, n); vl <- rep(0, n); mip <- rep(0, n)
  mpc <- rep(0, n); mlc <- rep(0, n); vlc <- rep(0, n); mipc <- rep(0, n)
  for (fit in fits) {
    nd <- ncol(fit$Wsamp); thin <- (fit$M - fit$burn) / nd
    didx <- fit$burn + seq_len(nd) * thin
    for (j in seq_len(nd)) {
      p <- fit$PAR[1:8, didx[j]]; names(p) <- rownames(fit$PAR)[1:8]
      lj <- log_marg(p, fit$PAR[9:12, didx[j]], f)
      lc <- lvc_logdens(q, y, fit$Wsamp[, j], p, f, 0.1)
      cnt <- cnt + 1
      mp <- mp + (exp(lj) - mp) / cnt; mip <- mip + (exp(-lj) - mip) / cnt
      d1 <- lj - ml_; ml_ <- ml_ + d1 / cnt; vl <- vl + d1 * (lj - ml_)
      mpc <- mpc + (exp(lc) - mpc) / cnt; mipc <- mipc + (exp(-lc) - mipc) / cnt
      d2 <- lc - mlc; mlc <- mlc + d2 / cnt; vlc <- vlc + d2 * (lc - mlc)
    }
  }
  lp <- sum(log(pmax(mp, 1e-300))); pw <- sum(vl / cnt)
  lpc <- sum(log(pmax(mpc, 1e-300))); pwc <- sum(vlc / cnt)
  crit <- rbind(crit, data.frame(family = f,
    lppd_m = lp, pwaic_m = pw, WAIC_m = -2 * lp + 2 * pw, CVML_m = -sum(log(mip)),
    WAIC_c = -2 * lpc + 2 * pwc, CVML_c = -sum(log(mipc))))
  cat(sprintf("%-13s marginal WAIC %8.1f  CVML %7.1f | conditional WAIC %9.1f\n",
              f, -2 * lp + 2 * pw, -sum(log(mip)), -2 * lpc + 2 * pwc))

  ## --- posterior summary ---------------------------------------------------
  S <- do.call(cbind, L$post)
  S2 <- S; S2["lalpha", ] <- exp(S["lalpha", ]); S2["lsigma", ] <- exp(S["lsigma", ])
  rownames(S2)[rownames(S2) == "lalpha"] <- "alpha"
  rownames(S2)[rownames(S2) == "lsigma"] <- "sigma"
  if (f == "Independence") S2 <- S2[!rownames(S2) %in% c("xi0", "xi1"), ]
  summ[[f]] <- t(apply(S2, 1, function(v)
    c(mean = mean(v), sd = sd(v), lo = quantile(v, .025, names = FALSE),
      hi = quantile(v, .975, names = FALSE))))
}

cat("\n=== convergence ===\n"); print(diagn, digits = 4, row.names = FALSE)

crit$dWAIC_m <- crit$WAIC_m - min(crit$WAIC_m)
crit$dWAIC_c <- crit$WAIC_c - min(crit$WAIC_c)
cat("\n=== model comparison ===\n")
print(crit[order(crit$WAIC_m), c("family", "WAIC_m", "dWAIC_m", "CVML_m",
                                 "WAIC_c", "dWAIC_c")], digits = 5, row.names = FALSE)
best <- crit$family[which.min(crit$WAIC_m)]
cat("\nselected by marginal WAIC:", best, "\n")
cat("selected by conditional WAIC:", crit$family[which.min(crit$WAIC_c)], "\n\n")
print(round(summ[[best]], 4))

## Both lambda's are positive and both markers rise as filtration falls, so W is
## oriented as renal impairment; P(xi1 > 0) measures the evidence that the
## dependence between the two markers strengthens with impairment.
for (f in setdiff(fams, "Independence")) {
  S <- do.call(cbind, load_fam(f)$post)
  cat(sprintf("P(xi1 > 0 | data), %-9s = %.4f\n", f, mean(S["xi1", ] > 0)))
}

## ------------------------------------------------------------- LaTeX tables -
esc <- function(x) formatC(x, format = "f", digits = 1)
tb <- crit[order(crit$WAIC_m), ]
writeLines(c(
  "\\begin{tabular}{lrrrrr}", "\\toprule",
  "& \\multicolumn{3}{c}{marginal} & \\multicolumn{2}{c}{conditional} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}",
  "copula & WAIC & $\\Delta$WAIC & CVML & WAIC & $\\Delta$WAIC \\\\", "\\midrule",
  paste0(tb$family, " & ", esc(tb$WAIC_m), " & ", esc(tb$dWAIC_m), " & ",
         esc(tb$CVML_m), " & ", esc(tb$WAIC_c), " & ", esc(tb$dWAIC_c), " \\\\"),
  "\\bottomrule", "\\end{tabular}"), "output/tables/model_comparison.tex")

sb <- summ[[best]]
writeLines(c(
  "\\begin{tabular}{lrrr}", "\\toprule",
  "parameter & mean & sd & 95\\% credible interval \\\\", "\\midrule",
  paste0(ifelse(rownames(sb) %in% names(pretty), pretty[rownames(sb)], rownames(sb)),
         " & ", formatC(sb[, "mean"], format = "f", digits = 3),
         " & ", formatC(sb[, "sd"], format = "f", digits = 3),
         " & $(", formatC(sb[, "lo"], format = "f", digits = 3), ",\\ ",
         formatC(sb[, "hi"], format = "f", digits = 3), ")$ \\\\"),
  "\\bottomrule", "\\end{tabular}"), "output/tables/posterior_best.tex")

## ------------------------------------------------------------------ figures -
base <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "plain"))

L <- load_fam(best); S <- do.call(cbind, L$post)
Wall <- do.call(cbind, lapply(L$fits, function(f) f$Wsamp))
wgrid <- seq(quantile(Wall, .01), quantile(Wall, .99), length.out = 120)
idx <- seq(1, ncol(S), length.out = 1200)
TH <- sapply(idx, function(j)
  theta_of_w(wgrid, c(S["xi0", j], S["xi1", j]), best, 0.1))
TA <- apply(TH, 2, tau_of_theta, family = best)

band <- function(M, lab, file) {
  d <- data.frame(w = wgrid, med = apply(M, 1, median),
                  lo = apply(M, 1, quantile, .025),
                  hi = apply(M, 1, quantile, .975))
  p <- ggplot(d, aes(w, med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey75", alpha = .6) +
    geom_line(linewidth = .7) +
    labs(x = expression("latent renal impairment "*W), y = lab) + base
  ggsave(file, p, width = 5.2, height = 3.8, dpi = 300)
}
band(TH, expression("copula parameter "*theta*"("*W*")"), "output/figs/NHANES_theta_of_W.png")
band(TA, expression("Kendall's "*tau*"("*W*")"), "output/figs/NHANES_tau_of_W.png")

Wm <- rowMeans(sapply(L$fits, function(f) f$Wmean))
dd <- data.frame(W = Wm, egfr = NH$dat$egfr, age = NH$dat$age)
p <- ggplot(dd, aes(egfr, W)) +
  geom_point(alpha = .12, size = .6) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "black", linewidth = .7) +
  labs(x = expression("CKD-EPI eGFR (mL/min/1.73"*m^2*")"),
       y = expression("posterior mean of "*W)) + base
ggsave("output/figs/NHANES_W_vs_eGFR.png", p, width = 5.2, height = 3.8, dpi = 300)
cat("\nSpearman(W, eGFR) =", round(cor(Wm, NH$dat$egfr, method = "spearman"), 3), "\n")

png("output/figs/NHANES_trace.png", width = 1900, height = 1300, res = 165)
op <- par(mfrow = c(3, 4), mar = c(3.4, 3.4, 2, 1), mgp = c(2.1, .7, 0))
for (v in rownames(S)) {
  plot(L$post[[1]][v, seq(1, ncol(L$post[[1]]), by = 50)], type = "l",
       col = "grey35", xlab = "iteration/50", ylab = "", main = v)
  lines(L$post[[2]][v, seq(1, ncol(L$post[[2]]), by = 50)], col = "steelblue")
}
par(op); dev.off()

saveRDS(list(crit = crit, summ = summ, diagn = diagn, best = best), "output/results.rds")
cat("\ntables/ and figs/ written; best model:", best, "\n")
