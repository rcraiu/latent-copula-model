## Is the faster proposal actually better? Speed per iteration is only useful if
## the chain does not mix correspondingly worse, so compare effective sample size
## per second rather than seconds per iteration.
source("R/lvc_fast.R")
suppressMessages(library(coda))

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X
M <- as.integer(Sys.getenv("BENCH_M", "20000"))

res <- NULL
for (fam in c("Gaussian", "Gumbel")) {
  for (bw in c(FALSE, TRUE)) {
    tm <- system.time(
      f <- run_lvc_fast(q, y, X, family = fam, M = M, burn = M / 2, seed = 21,
                        verbose = 0, track_criteria = FALSE, blockwise = bw)
    )[["elapsed"]]
    S <- f$PAR[, (f$burn + 1):f$M]
    keep <- if (fam == "Independence") setdiff(rownames(S), c("xi0", "xi1")) else rownames(S)
    e <- effectiveSize(mcmc(t(S[keep, ])))
    res <- rbind(res, data.frame(
      family = fam, proposal = if (bw) "joint block" else "componentwise",
      sec = tm, minESS = min(e), medESS = median(e),
      minESS_per_sec = min(e) / tm))
    cat(sprintf("%-9s %-14s %6.1fs  min ESS %7.1f  median ESS %8.1f  min ESS/sec %6.2f\n",
                fam, if (bw) "joint block" else "componentwise", tm, min(e),
                median(e), min(e) / tm))
    cat("    worst-mixing parameters: ",
        paste(names(sort(e))[1:3], round(sort(e)[1:3], 1), collapse = ", "), "\n")
  }
}
cat("\n")
for (fam in unique(res$family)) {
  r <- res[res$family == fam, ]
  cat(sprintf("%s: joint block gives %.2fx the min ESS per second\n", fam,
              r$minESS_per_sec[2] / r$minESS_per_sec[1]))
}
saveRDS(res, "output/bench_ess.rds")
