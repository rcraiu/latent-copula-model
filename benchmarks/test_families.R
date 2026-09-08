## Validation of the copula families, in two parts:
##   1. the per-observation log density, R against C++, at random parameter
##      states, which must agree to machine precision;
##   2. the full sampler path, R against C++, which must agree because the two
##      consume the random number stream in the same order.
##
## Part 1 is the more sensitive of the two. The samplers are chaotic: a
## difference of one bit in a density can flip a Metropolis decision and send
## the paths apart, so part 2 failing while part 1 passes indicates a numerical
## discrepancy rather than a wrong formula. This is exactly how the gamma
## rate/scale inconsistency was found -- R parameterised the gamma margin by
## rate and C++ by scale, and since 1/(alpha/b) and b/alpha differ in the last
## bit, the tails of pgamma disagreed by enough to matter.

suppressMessages(library(Rcpp))
source("R/lvc_gamma_gauss.R")
source("R/lvc_fast.R")
sourceCpp("R/lvc_kernel.cpp")

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X

FAMS <- c("Gaussian", "Clayton", "Frank", "Gumbel", "Joe", "SurvClayton")

## Tolerance for part 1. Every family should agree to machine precision except
## the Gaussian, where the two implementations are deliberately not identical:
## the normal margin makes Phi^{-1}(F_y(y|w)) the identity, which the C++ code
## uses directly while the R reference evaluates qnorm(pnorm(z)) and loses
## precision in the tails. The C++ value is the more accurate one, so this is a
## known and accepted discrepancy rather than a defect.
DENS_TOL <- c(Gaussian = 1e-6, Clayton = 1e-11, Frank = 1e-11,
              Gumbel = 1e-11, Joe = 1e-11, SurvClayton = 1e-11)

## -- 1. densities -------------------------------------------------------
cat("--- per-observation log density, R vs C++, 30 random states each ---\n")
set.seed(3)
ok_all <- TRUE
for (fam in FAMS) {
  worst <- 0
  for (rep in 1:30) {
    W <- rnorm(length(q))
    p <- c(muq = rnorm(1, -.05, .1), lamq = runif(1, .05, .3),
           lalpha = log(runif(1, 20, 300)),
           muy = rnorm(1, -.2, .1), lamy = runif(1, .05, .3),
           lsigma = log(runif(1, .2, .4)),
           xi0 = rnorm(1, 0, 4), xi1 = rnorm(1, 1, 1))
    rr <- lvc_logdens(q, y, W, p, fam, 0.1, "joint")
    cc <- lvc_logdens_cpp(q, y, W, p[["muq"]], p[["lamq"]], exp(p[["lalpha"]]),
                          p[["muy"]], p[["lamy"]], exp(p[["lsigma"]]),
                          p[["xi0"]], p[["xi1"]], FAMILY_CODE[[fam]], 0.1, 0L)
    keep <- is.finite(rr) & is.finite(cc)
    worst <- max(worst, max(abs((rr[keep] - cc[keep]) /
                                pmax(abs(rr[keep]), 1))))
  }
  pass <- worst < DENS_TOL[[fam]]
  ok_all <- ok_all && pass
  cat(sprintf("%-12s max relative difference %.3e  %s%s\n", fam, worst,
              if (pass) "ok" else "*** FAIL ***",
              if (fam == "Gaussian") "  (qnorm(pnorm(z)) in R; see comment)" else ""))
}

## -- 2. sampler paths ---------------------------------------------------
cat("\n--- full sampler path, R vs C++, componentwise, M = 1500 ---\n")
M <- 1500
for (fam in FAMS) {
  a <- run_lvc     (q, y, X, family = fam, M = M, burn = M / 2, seed = 7,
                    verbose = 0, track_criteria = FALSE)
  b <- run_lvc_fast(q, y, X, family = fam, M = M, burn = M / 2, seed = 7,
                    verbose = 0, track_criteria = FALSE, blockwise = FALSE)
  d <- max(abs(a$PAR - b$PAR))
  pass <- d < 1e-6
  ok_all <- ok_all && pass
  cat(sprintf("%-12s max absolute difference %.3e  %s\n", fam, d,
              if (pass) "ok" else "*** FAIL ***"))
}

## -- 3. against VineCopula ---------------------------------------------
if (requireNamespace("VineCopula", quietly = TRUE)) {
  cat("\n--- copula densities against VineCopula at fixed parameters ---\n")
  set.seed(1)
  u1 <- runif(200, .01, .99); u2 <- runif(200, .01, .99)
  ## VineCopula codes: 1 Gaussian, 3 Clayton, 4 Gumbel, 5 Frank, 6 Joe,
  ## 13 survival (180 degree) Clayton
  spec <- list(list("Gaussian", 1, .6), list("Clayton", 3, 2), list("Gumbel", 4, 1.8),
               list("Frank", 5, 4), list("Joe", 6, 2), list("SurvClayton", 13, 2))
  for (s in spec) {
    mine <- log_copula_density(u1, u2, rep(s[[3]], 200), s[[1]])
    theirs <- log(VineCopula::BiCopPDF(u1, u2, family = s[[2]], par = s[[3]]))
    d <- max(abs(mine - theirs))
    ok_all <- ok_all && d < 1e-10
    cat(sprintf("%-12s par %4.1f  max absolute difference %.3e  %s\n",
                s[[1]], s[[3]], d, if (d < 1e-10) "ok" else "*** FAIL ***"))
  }
}

cat("\n", if (ok_all) "all checks passed" else "*** SOME CHECKS FAILED ***", "\n")
