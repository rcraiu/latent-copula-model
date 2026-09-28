## Compare the constrained (||beta|| = 1) and unconstrained Frank fits on the
## NHANES kidney data: scale-ridge correlations in the free fit, and effective
## sample sizes under the two parametrisations at identical run length.
##
## Expects output/fits/Frank_chain{1,2}.rds   (beta = norm1)
##         output/fits/Frank_chain{11,12}.rds (beta = free)

suppressMessages(library(coda))

load_pair <- function(chains) {
  lapply(chains, function(k) {
    fit <- readRDS(sprintf("output/fits/Frank_chain%d.rds", k))
    fit$PAR[, (fit$burn + 1):fit$M]
  })
}

norm1 <- load_pair(c(1, 2))
free  <- load_pair(c(11, 12))

## ---- scale ridge in the unconstrained fit ----------------------------------
Sf <- do.call(cbind, free)
bnorm <- sqrt(colSums(Sf[c("beta1", "beta2", "beta3", "beta4"), ]^2))
cat("=== unconstrained Frank fit (chains 11, 12 pooled) ===\n")
cat(sprintf("||beta||: mean %.3f  sd %.3f  95%% CI (%.3f, %.3f)\n",
            mean(bnorm), sd(bnorm),
            quantile(bnorm, .025), quantile(bnorm, .975)))
for (v in c("lamq", "lamy", "xi1", "xi0", "lalpha"))
  cat(sprintf("cor(log||beta||, %-6s) = %6.3f\n", v, cor(log(bnorm), Sf[v, ])))

## ---- ESS comparison ---------------------------------------------------------
ess_of <- function(post) {
  ml <- mcmc.list(lapply(post, function(S)
    mcmc(t(S[, seq(1, ncol(S), by = 10)]))))
  effectiveSize(ml)
}
e1 <- ess_of(norm1); e0 <- ess_of(free)
cat("\n=== ESS, two chains pooled (thin 10), same run length ===\n")
tab <- data.frame(free = round(e0[names(e1)]), norm1 = round(e1),
                  ratio = round(e1 / e0[names(e1)], 2))
print(tab)
cat(sprintf("\nmin ESS: free %d   norm1 %d\n",
            round(min(e0)), round(min(e1))))

## ---- posterior agreement on identified quantities ---------------------------
Sn <- do.call(cbind, norm1)
cat("\n=== identified quantities: free vs norm1 ===\n")
u_free <- Sf[c("beta1","beta2","beta3","beta4"), ] /
          rep(bnorm, each = 4)                      # direction in free fit
for (j in 1:4)
  cat(sprintf("beta%d direction: free %6.3f   norm1 %6.3f\n",
              j, mean(u_free[j, ]), mean(Sn[paste0("beta", j), ])))
for (v in c("lamq", "lamy", "xi1"))
  cat(sprintf("%s: free %.4f x mean||beta|| %.3f = %.4f   norm1 %.4f\n",
              v, mean(Sf[v, ]), mean(bnorm), mean(Sf[v, ] * bnorm),
              mean(Sn[v, ])))
