## Validate run_lvc_fast() against run_lvc() and measure the speedup.
source("R/lvc_gamma_gauss.R")     # pure R reference
source("R/lvc_fast.R")            # Rcpp driver

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X

M <- as.integer(Sys.getenv("TEST_M", "2000"))
cat("n =", length(q), "  M =", M, "\n\n")

for (fam in c("Gaussian", "Clayton", "Frank", "Gumbel", "Independence")) {
  tR <- system.time(
    a <- run_lvc(q, y, X, family = fam, M = M, burn = M / 2, seed = 99,
                 verbose = 0, thin_W = 0))[["elapsed"]]
  tC <- system.time(
    b <- run_lvc_fast(q, y, X, family = fam, M = M, burn = M / 2, seed = 99,
                      verbose = 0, thin_W = 0))[["elapsed"]]

  dpar <- max(abs(a$PAR - b$PAR))
  dW   <- max(abs(a$W - b$W))
  dcr  <- if (!is.null(a$criteria)) max(abs(a$criteria - b$criteria)) else NA

  cat(sprintf("%-13s R %6.2fs  C++ %6.2fs  (%5.1fx)  | max|dPAR| %.2e  max|dW| %.2e  max|dcrit| %.2e\n",
              fam, tR, tC, tR / tC, dpar, dW, dcr))
}

cat("\n--- joint block proposals instead of componentwise ---\n")
for (fam in c("Gaussian", "Gumbel")) {
  tb <- system.time(
    bb <- run_lvc_fast(q, y, X, family = fam, M = M, burn = M / 2, seed = 99,
                       verbose = 0, blockwise = TRUE))[["elapsed"]]
  cat(sprintf("%-13s blockwise %6.2fs   acceptance: %s\n", fam, tb,
              paste(sprintf("%.2f", bb$acc_blocks[1:6]), collapse = " ")))
}
