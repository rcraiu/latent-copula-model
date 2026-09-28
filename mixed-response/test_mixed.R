## Checks for mixed_model.R. Run from the repository root:
##   Rscript mixed-response/test_mixed.R
##
## 1. h-functions against VineCopula at fixed parameter values
## 2. hinv(h(u1 | u2) | u2) recovers u1 for every family
## 3. the Gauss-Hermite marginalisation of W against the closed-form margin of
##    y and against a brute-force fine grid for the joint density

source("mixed-response/mixed_model.R")
suppressMessages(library(VineCopula))
set.seed(1)

n  <- 400
u1 <- runif(n, 0.005, 0.995); u2 <- runif(n, 0.005, 0.995)
vf <- c(Gaussian = 1, Clayton = 3, Gumbel = 4, Frank = 5)
pars <- list(Gaussian = c(-0.7, 0.3, 0.9), Clayton = c(0.3, 2, 8),
             Gumbel = c(1.2, 2.5, 6), Frank = c(-8, 0.5, 12))

cat("--- h-functions against VineCopula ---\n")
worst <- 0
for (f in names(vf)) for (par in pars[[f]]) {
  mine   <- hfunc(u1, u2, rep(par, n), f)
  theirs <- BiCopHfunc2(u1, u2, vf[[f]], par)
  e <- max(abs(mine - theirs)); worst <- max(worst, e)
  cat(sprintf("  %-9s theta = %5.2f  max abs diff %.2e\n", f, par, e))
}
stopifnot(worst < 1e-8)

cat("--- inverse h-functions (round trip), varying theta ---\n")
for (f in c(names(vf), "Independence")) {
  th <- switch(f, Gaussian = tanh(rnorm(n)), Clayton = exp(rnorm(n)),
               Gumbel = 1 + exp(rnorm(n)), Frank = 3 * rnorm(n), Independence = rep(0, n))
  v  <- hfunc(u1, u2, th, f)
  ## where h saturates at 0 or 1 in double precision the inverse is undefined
  ok <- v > 1e-9 & v < 1 - 1e-9
  e  <- max(abs(hinv(v, u2, th, f) - u1)[ok])
  cat(sprintf("  %-13s max abs error %.2e  (%d of %d non-saturated)\n", f, e, sum(ok), n))
  stopifnot(e < 1e-6)
}

cat("--- large-parameter stability ---\n")
for (f in names(vf)) {
  th <- switch(f, Gaussian = rep(0.999999, n), Clayton = rep(150, n),
               Gumbel = rep(150, n), Frank = rep(-150, n))
  v <- hfunc(u1, u2, th, f)
  cat(sprintf("  %-9s finite: %s  in [0,1]: %s\n", f, all(is.finite(v)), all(v >= 0 & v <= 1)))
  stopifnot(all(is.finite(v)))
}

cat("--- W marginalisation ---\n")
p  <- c(mud = -1, lamd = 1, muc = 0, lamc = 1, ssq = 1)
xi <- c(-1, 3)
X  <- cbind(runif(n, -3, 1), runif(n, -1, 1)); X <- scale(X)
m  <- as.numeric(X %*% c(1, 0.5))
W  <- m + rnorm(n)
for (f in c("Gaussian", "Gumbel", "Clayton", "Frank", "Independence")) {
  d  <- simulate_mixed(W, p, xi, f)
  L  <- logdens_marg_W(d$q, d$y, m, p, xi, f)
  ## the same integrals on a grid ten times finer
  R  <- logdens_marg_W(d$q, d$y, m, p, xi, f, gh = simpson_nodes(8, 0.002))
  ## the margin of y implied by the quadrature, joint / p(q | y), against its
  ## closed form; this checks the joint and the ratio together
  ly_quad <- L[, "joint"] - L[, "d_g_c"]
  ly_true <- dnorm(d$y, p["muc"] + p["lamc"] * m, sqrt(p["ssq"] + p["lamc"]^2), log = TRUE)
  ## the margin of q from the quadrature against the same grid refinement
  cat(sprintf("  %-13s  max |err| joint %.1e  p(q|y) %.1e  p(y|q) %.1e   closed-form p(y) %.1e   max q %d\n",
              f, max(abs(L[, 1] - R[, 1])), max(abs(L[, 2] - R[, 2])),
              max(abs(L[, 3] - R[, 3])), max(abs(ly_quad - ly_true)), max(d$q)))
  stopifnot(max(abs(L - R)) < 1e-5, max(abs(ly_quad - ly_true)) < 1e-10)
  ## Kendall's tau of the simulated pair, as a sanity check on the direction
  cat(sprintf("                 Kendall tau(q, y) = %.3f, theta range %.2f to %.2f\n",
              cor(d$q, d$y, method = "kendall"), min(d$theta), max(d$theta)))
}
cat("\nall checks passed\n")
