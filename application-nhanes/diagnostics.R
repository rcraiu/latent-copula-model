## Convergence diagnostics and posterior summaries across the two chains.
suppressMessages(library(coda))

fams <- c("Gaussian", "Clayton", "Frank", "Gumbel", "Independence")
tr <- function(p) {                       # report on the interpretable scale
  out <- p
  rownames(out)[rownames(out) == "lalpha"] <- "alpha"
  rownames(out)[rownames(out) == "lsigma"] <- "sigma"
  out["alpha", ] <- exp(p["lalpha", ]); out["sigma", ] <- exp(p["lsigma", ])
  out
}

for (f in fams) {
  ch <- lapply(1:2, function(k) {
    fit <- readRDS(sprintf("output/fits/%s_chain%d.rds", f, k))
    S <- tr(fit$PAR[, (fit$burn + 1):fit$M])
    if (f == "Independence") S <- S[!rownames(S) %in% c("xi0", "xi1"), ]
    mcmc(t(S[, seq(1, ncol(S), by = 10)]))
  })
  ml <- mcmc.list(ch)
  g <- gelman.diag(ml, multivariate = TRUE, autoburnin = FALSE)
  ess <- effectiveSize(ml)
  cat(sprintf("\n===== %s =====\n", f))
  cat("max PSRF:", round(max(g$psrf[, 1]), 4),
      "  MPSRF:", round(g$mpsrf, 4),
      "  min ESS:", round(min(ess)), "\n")
  bad <- g$psrf[g$psrf[, 1] > 1.05, , drop = FALSE]
  if (nrow(bad)) { cat("  PSRF > 1.05 for:\n"); print(round(bad, 3)) }

  S <- do.call(cbind, lapply(ch, function(m) t(as.matrix(m))))
  su <- t(apply(S, 1, function(v)
    c(mean = mean(v), sd = sd(v), q025 = quantile(v, .025, names = FALSE),
      q500 = quantile(v, .5, names = FALSE), q975 = quantile(v, .975, names = FALSE))))
  print(round(su, 4))
  if (f != "Independence")
    cat("P(xi1 > 0 | data) =", round(mean(S["xi1", ] > 0), 4), "\n")
}
