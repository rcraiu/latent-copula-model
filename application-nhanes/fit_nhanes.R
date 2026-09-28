## Fit the gamma-Gaussian latent variable copula model to the NHANES kidney
## data, one copula family and one chain per invocation:
##     Rscript fit_nhanes.R <family> <chain> [M] [engine] [proposal] [beta]
##
## engine   "cpp" (default) uses the Rcpp sampler in lvc_sampler.cpp;
##          "r" uses the reference implementation in lvc_gamma_gauss.R.
## proposal "block" (default) proposes each parameter block jointly from the
##          adapted covariance and appends a joint (alpha, xi) update;
##          "blocknox" drops that extra update; "component" updates one
##          coordinate at a time as in the reference implementation.
## beta     "norm1" (default) restricts beta to the unit sphere, frees the
##          residual variance of W, and uses the PX-DA group move of Xu et al.
##          (2016, JCGS, eq. 3.2); "free" keeps the unconstrained conjugate
##          draw with var(W|x) = 1.
##
## Measured min ESS per second relative to "component" at M = 20000:
##                        joint blocks    + joint (alpha, xi)
##          Gaussian          2.73x              3.58x
##          Gumbel            1.95x              2.73x
##          Frank             3.97x              1.81x
## (the Frank ordering is within the noise of an ESS estimate this small).

args     <- commandArgs(trailingOnly = TRUE)
family   <- args[1]
chain    <- as.integer(args[2])
M        <- if (length(args) >= 3) as.integer(args[3]) else 60000
engine   <- if (length(args) >= 4) args[4] else "cpp"
proposal <- if (length(args) >= 5) args[5] else "block"
betaopt  <- if (length(args) >= 6) args[6] else "norm1"
bnorm1   <- (betaopt == "norm1")

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q  <- NH$q          # 1 / cystatin C               -> gamma margin
y  <- NH$y          # 1 / beta-2 microglobulin     -> normal margin
X  <- NH$X          # age, female, BMI, SBP (standardised, no intercept)

dir.create("fits", showWarnings = FALSE)
out <- sprintf("output/fits/%s_chain%d.rds", family, chain)

cat(sprintf("fitting %s, chain %d, n = %d, M = %d, engine = %s, proposal = %s, beta = %s\n",
            family, chain, length(q), M, engine, proposal, betaopt))

if (engine == "cpp") {
  source("R/lvc_fast.R")
  tm <- system.time(
    fit <- run_lvc_fast(q, y, X, family = family, M = M, burn = floor(M / 2),
                        gscale = 0.1, seed = 1729 + 1000 * chain,
                        verbose = M / 6, thin_W = floor(M / 2 / 300),
                        blockwise = (proposal %in% c("block", "blocknox")),
                        cross_block = (proposal == "block"),
                        beta_norm1 = bnorm1, wvar_free = bnorm1, px_da = bnorm1)
  )
} else {
  source("R/lvc_gamma_gauss.R")
  tm <- system.time(
    fit <- run_lvc(q, y, X, family = family, M = M, burn = floor(M / 2),
                   gscale = 0.1, seed = 1729 + 1000 * chain,
                   verbose = M / 6, thin_W = floor(M / 2 / 300),
                   beta_norm1 = bnorm1, wvar_free = bnorm1, px_da = bnorm1)
  )
}
fit$elapsed <- tm[["elapsed"]]
fit$chain <- chain
fit$engine <- engine
fit$proposal <- proposal
fit$beta_constraint <- betaopt
saveRDS(fit, out)

cat(sprintf("\ndone in %.1f min -> %s\n", tm[["elapsed"]] / 60, out))
if (!is.null(fit$criteria)) print(round(fit$criteria, 1))
cat("\nacceptance rates:\n")
if (engine == "cpp") {
  print(round(fit$acc_blocks, 3)); cat("W:", round(fit$acc_W, 3), "\n")
} else {
  print(lapply(fit$acc$blocks, round, 3)); cat("W:", round(fit$acc$W, 3), "\n")
}

S <- fit$PAR[, (fit$burn + 1):fit$M, drop = FALSE]
if (family == "Independence") S <- S[!rownames(S) %in% c("xi0", "xi1"), ]
S <- rbind(S, alpha = exp(S["lalpha", ]), sigma2 = exp(2 * S["lsigma", ]))
print(round(t(apply(S, 1, function(v)
  c(mean = mean(v), sd = sd(v), quantile(v, c(0.025, 0.5, 0.975))))), 4))
