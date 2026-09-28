## Validate the beta_norm1 (unit-sphere beta) path of run_lvc_fast() against
## the pure R reference, on the NHANES data. Run from the repository root.
source("R/lvc_gamma_gauss.R")
source("R/lvc_fast.R")

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X

M <- as.integer(Sys.getenv("TEST_M", "2000"))
cat("n =", length(q), "  M =", M, "\n\n")

for (fam in c("Gaussian", "Gumbel", "Independence")) {
  tR <- system.time(
    a <- run_lvc(q, y, X, family = fam, M = M, burn = M / 2, seed = 99,
                 verbose = 0, thin_W = 0, beta_norm1 = TRUE))[["elapsed"]]
  tC <- system.time(
    b <- run_lvc_fast(q, y, X, family = fam, M = M, burn = M / 2, seed = 99,
                      verbose = 0, thin_W = 0, beta_norm1 = TRUE))[["elapsed"]]

  dpar <- max(abs(a$PAR - b$PAR))
  dW   <- max(abs(a$W - b$W))
  bn   <- sqrt(colSums(b$PAR[9:nrow(b$PAR), ]^2))
  cat(sprintf("%-13s R %6.2fs  C++ %6.2fs (%5.1fx) | max|dPAR| %.2e  max|dW| %.2e | ||beta|| in [%.6f, %.6f]  acc_beta R %.3f C++ %.3f\n",
              fam, tR, tC, tR / tC, dpar, dW, min(bn), max(bn),
              a$acc$beta, b$acc_beta))
}
