## Does the joint (alpha, xi) block actually fix the mixing problem?
## Compare effective sample size per second across the three configurations.
source("R/lvc_fast.R")
suppressMessages(library(coda))

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X
M <- as.integer(Sys.getenv("BENCH_M", "20000"))

cfg <- list(
  list(lab = "componentwise",           bw = FALSE, cb = FALSE),
  list(lab = "joint blocks",            bw = TRUE,  cb = FALSE),
  list(lab = "joint blocks + alpha-xi", bw = TRUE,  cb = TRUE)
)

res <- NULL
for (fam in c("Gaussian", "Gumbel", "Frank")) {
  cat("\n===", fam, "===\n")
  for (cf in cfg) {
    tm <- system.time(
      f <- run_lvc_fast(q, y, X, family = fam, M = M, burn = M / 2, seed = 21,
                        verbose = 0, track_criteria = FALSE,
                        blockwise = cf$bw, cross_block = cf$cb)
    )[["elapsed"]]
    S <- f$PAR[, (f$burn + 1):f$M]
    e <- effectiveSize(mcmc(t(S)))
    res <- rbind(res, data.frame(family = fam, config = cf$lab, sec = tm,
                                 minESS = min(e), ESS_alpha = e[["lalpha"]],
                                 ESS_xi0 = e[["xi0"]],
                                 minESS_per_sec = min(e) / tm))
    cat(sprintf("  %-24s %6.1fs  minESS %6.1f  ESS(alpha) %6.1f  ESS(xi0) %6.1f  minESS/s %5.3f\n",
                cf$lab, tm, min(e), e[["lalpha"]], e[["xi0"]], min(e) / tm))
    cat("     slowest: ", paste(names(sort(e))[1:3], round(sort(e)[1:3], 1),
                                collapse = ", "), "\n")
    ## posterior correlation between alpha and the copula parameters
    if (cf$cb) cat(sprintf("     cor(log alpha, xi0) = %.3f\n",
                           cor(S["lalpha", ], S["xi0", ])))
  }
}

cat("\n=== gain over componentwise, in min ESS per second ===\n")
for (fam in unique(res$family)) {
  r <- res[res$family == fam, ]
  cat(sprintf("%-9s joint blocks %.2fx    joint blocks + alpha-xi %.2fx\n", fam,
              r$minESS_per_sec[2] / r$minESS_per_sec[1],
              r$minESS_per_sec[3] / r$minESS_per_sec[1]))
}
saveRDS(res, "output/bench_crossblock.rds")
