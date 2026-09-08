## Parameter-recovery check for the gamma-Gaussian latent variable copula
## sampler: simulate from the model of Section 2.1 and see whether the
## posterior concentrates near the truth.

source("R/lvc_gamma_gauss.R")

set.seed(20260907)

n  <- 2000
p  <- 4
gs <- 0.1

beta_true <- c(0.6, -0.4, 0.3, 0.5); beta_true <- beta_true / sqrt(sum(beta_true^2))
muq_true  <- 0.00; lamq_true <- 0.40; alpha_true <- 20
muy_true  <- -0.20; lamy_true <- 0.30; sig_true  <- 0.25
## kept mild so that tanh(0.1*(xi0+xi1*w)) does not saturate at +/-1
xi_true   <- c(2, 4)

X <- matrix(rnorm(n * p), n, p)
X <- scale(X); X <- matrix(as.numeric(X), n, p)
W <- as.numeric(X %*% beta_true) + rnorm(n)

rho <- theta_of_w(W, xi_true, "Gaussian", gs)
z1  <- rnorm(n)
z2  <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
u1  <- pnorm(z1); u2 <- pnorm(z2)

mq <- exp(muq_true + lamq_true * W)
q  <- qgamma(u1, shape = alpha_true, rate = alpha_true / mq)
y  <- qnorm(u2, mean = muy_true + lamy_true * W, sd = sig_true)

cat("simulated data: n =", n, "\n")
cat(sprintf("  q: median %.2f  range [%.2f, %.2f]\n", median(q), min(q), max(q)))
cat(sprintf("  y: median %.2f  range [%.2f, %.2f]\n", median(y), min(y), max(y)))
cat("  implied rho range: [", paste(round(range(rho), 3), collapse = ", "), "]\n")
cat("  Kendall tau =", round(cor(q, y, method = "kendall"), 3), "\n\n")

M <- as.integer(Sys.getenv("LVC_M", "40000"))
tm <- system.time(
  fit <- run_lvc(q, y, X, family = "Gaussian", M = M, burn = floor(M / 2),
                 gscale = gs, seed = 7, verbose = M / 10, track_criteria = FALSE)
)
cat("\nelapsed:", round(tm[["elapsed"]], 1), "sec for", M, "iterations\n\n")

S <- lvc_summary(fit)
truth <- c(muq = muq_true, lamq = lamq_true, lalpha = log(alpha_true),
           muy = muy_true, lamy = lamy_true, lsigma = log(sig_true),
           xi0 = xi_true[1], xi1 = xi_true[2],
           beta1 = beta_true[1], beta2 = beta_true[2],
           beta3 = beta_true[3], beta4 = beta_true[4],
           alpha = alpha_true, sigma2 = sig_true^2)
S <- cbind(truth = truth[rownames(S)], round(S, 4))
print(S)

covered <- truth[rownames(S)] >= S[, "2.5%"] & truth[rownames(S)] <= S[, "97.5%"]
cat("\n95% credible interval covers the truth for",
    sum(covered), "of", length(covered), "parameters\n")
print(names(covered)[!covered])

cat("\nacceptance rates:\n")
print(lapply(fit$acc$blocks, round, 3))
cat("W:", round(fit$acc$W, 3), "\n")
