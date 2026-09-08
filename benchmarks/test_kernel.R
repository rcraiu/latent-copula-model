## Validate the Rcpp kernel against the pure-R implementation and time both.
suppressMessages(library(Rcpp))
source("R/copula_densities.R")
source("R/lvc_gamma_gauss.R")
sourceCpp("R/lvc_kernel.cpp")

FAM <- c(Gaussian = 0L, Clayton = 1L, Frank = 2L, Gumbel = 3L, Independence = 4L)

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; n <- length(q)
set.seed(4)
W <- rnorm(n)

p <- c(muq = -0.05, lamq = 0.2, lalpha = log(90), muy = -0.2, lamy = 0.15,
       lsigma = log(0.3), xi0 = 2.5, xi1 = 1.2)

cpp_call <- function(fam, W, p, what = 0L)
  lvc_logdens_cpp(q, y, W, p[["muq"]], p[["lamq"]], exp(p[["lalpha"]]),
                  p[["muy"]], p[["lamy"]], exp(p[["lsigma"]]),
                  p[["xi0"]], p[["xi1"]], FAM[[fam]], 0.1, what)

cat("=== agreement with the pure-R implementation (n =", n, ") ===\n")
for (fam in names(FAM)) {
  for (w_ in 0:2) {
    what <- c("joint", "q_g_y", "y_g_q")[w_ + 1]
    a <- lvc_logdens(q, y, W, p, fam, 0.1, what)
    b <- cpp_call(fam, W, p, w_)
    fin <- is.finite(a) & is.finite(b)
    cat(sprintf("%-13s %-6s max abs diff %.3e   (finite %d/%d, agree on -Inf: %s)\n",
                fam, what, max(abs(a[fin] - b[fin])), sum(fin), n,
                identical(is.finite(a), is.finite(b))))
  }
}

cat("\n=== the three-in-one call matches the three separate calls ===\n")
for (fam in names(FAM)) {
  M3 <- lvc_logdens3_cpp(q, y, W, p[["muq"]], p[["lamq"]], exp(p[["lalpha"]]),
                         p[["muy"]], p[["lamy"]], exp(p[["lsigma"]]),
                         p[["xi0"]], p[["xi1"]], FAM[[fam]], 0.1)
  d <- max(abs(M3[, 1] - cpp_call(fam, W, p, 0L)),
           abs(M3[, 2] - cpp_call(fam, W, p, 1L)),
           abs(M3[, 3] - cpp_call(fam, W, p, 2L)))
  cat(sprintf("%-13s max abs diff %.3e\n", fam, d))
}

cat("\n=== loglik sum ===\n")
for (fam in names(FAM)) {
  a <- sum(lvc_logdens(q, y, W, p, fam, 0.1))
  b <- lvc_loglik_cpp(q, y, W, p[["muq"]], p[["lamq"]], exp(p[["lalpha"]]),
                      p[["muy"]], p[["lamy"]], exp(p[["lsigma"]]),
                      p[["xi0"]], p[["xi1"]], FAM[[fam]], 0.1)
  cat(sprintf("%-13s R %14.6f   C++ %14.6f   diff %.3e\n", fam, a, b, abs(a - b)))
}

## --- extreme values, where the clamping matters ---------------------------
cat("\n=== behaviour in the tails ===\n")
Wx <- c(-8, -5, 0, 5, 8, 12)
qx <- c(1e-6, 0.05, 0.9, 3, 12, 40); yx <- c(-4, -2, 0, 1.5, 4, 9)
for (fam in names(FAM)) {
  a <- lvc_logdens(qx, yx, Wx, p, fam, 0.1)
  b <- lvc_logdens_cpp(qx, yx, Wx, p[["muq"]], p[["lamq"]], exp(p[["lalpha"]]),
                       p[["muy"]], p[["lamy"]], exp(p[["lsigma"]]),
                       p[["xi0"]], p[["xi1"]], FAM[[fam]], 0.1, 0L)
  fin <- is.finite(a) & is.finite(b)
  cat(sprintf("%-13s max abs diff %.3e over %d finite of %d\n",
              fam, if (any(fin)) max(abs(a[fin] - b[fin])) else 0, sum(fin), length(qx)))
}

## --- timing ---------------------------------------------------------------
cat("\n=== timing, 300 evaluations at n =", n, "===\n")
for (fam in c("Gaussian", "Gumbel", "Independence")) {
  tr <- system.time(for (i in 1:300) lvc_logdens(q, y, W, p, fam, 0.1))[["elapsed"]]
  tc <- system.time(for (i in 1:300) cpp_call(fam, W, p, 0L))[["elapsed"]]
  ts <- system.time(for (i in 1:300)
    lvc_loglik_cpp(q, y, W, p[["muq"]], p[["lamq"]], exp(p[["lalpha"]]),
                   p[["muy"]], p[["lamy"]], exp(p[["lsigma"]]),
                   p[["xi0"]], p[["xi1"]], FAM[[fam]], 0.1))[["elapsed"]]
  cat(sprintf("%-13s R %6.3fs   C++ vector %6.3fs (%4.1fx)   C++ sum %6.3fs (%4.1fx)\n",
              fam, tr, tc, tr / tc, ts, tr / ts))
}
