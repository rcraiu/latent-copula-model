## Two questions:
##  (1) Is a gamma margin defensible for cystatin C *conditional on the latent
##      level*? The covariate-only GLM is a pessimistic proxy, because the
##      latent variable W absorbs most of the marker variation. We mimic
##      conditioning on W by conditioning on the OTHER marker plus covariates.
##  (2) Choose a clinically principled exclusion rather than a marker cutoff
##      picked to improve fit: drop kidney-failure (CKD stage 5) subjects.

suppressMessages({ library(dplyr) })

NH  <- readRDS("application-nhanes/nhanes_kidney.rds")
dat <- NH$dat
X   <- NH$X

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3

## --- CKD-EPI 2021 creatinine equation (race-free), for staging only ------
egfr_cr <- function(scr, age, female) {
  kap <- ifelse(female, 0.7, 0.9)
  alp <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / kap, 1)^alp * pmax(scr / kap, 1)^(-1.200) *
    0.9938^age * ifelse(female, 1.012, 1)
}
dat$egfr <- egfr_cr(dat$creat, dat$age, dat$female == 1)

cat("--- CKD stage distribution (CKD-EPI 2021 creatinine) ---\n")
stage <- cut(dat$egfr, c(-Inf, 15, 30, 45, 60, 90, Inf),
             labels = c("G5 (<15)", "G4 (15-30)", "G3b (30-45)",
                        "G3a (45-60)", "G2 (60-90)", "G1 (>=90)"))
print(table(stage))

keep <- dat$egfr >= 15
cat("\nexcluding kidney failure (eGFR < 15):", sum(!keep), "dropped, n =", sum(keep), "\n")

## --- (1) gamma vs lognormal, conditioning on a proxy for W ---------------
fit_families <- function(v, Z) {
  gl  <- glm(v ~ Z, family = Gamma(link = "log"))
  mu  <- fitted(gl)
  nll <- function(la) -sum(dgamma(v, shape = exp(la), rate = exp(la) / mu, log = TRUE))
  la  <- optimize(nll, c(-3, 8))$minimum
  k   <- ncol(as.matrix(Z)) + 2
  gA  <- 2 * nll(la) + 2 * k

  lm1 <- lm(log(v) ~ Z)
  s   <- sqrt(mean(resid(lm1)^2))
  ll  <- sum(dnorm(log(v), mean = fitted(lm1), sd = s, log = TRUE)) - sum(log(v))
  lA  <- -2 * ll + 2 * k

  c(gamma.AIC = gA, alpha = exp(la), lognormal.AIC = lA,
    diff.per.obs = (gA - lA) / length(v))
}

Xk <- X[keep, ]; dk <- dat[keep, ]

cat("\n=== cystatin C: gamma vs lognormal, by conditioning set (n =", nrow(dk), ") ===\n")
cat("\ncovariates only (age, sex, BMI, SBP):\n")
print(round(fit_families(dk$q, Xk), 3))
cat("\ncovariates + log creatinine (a proxy for the latent level):\n")
print(round(fit_families(dk$q, cbind(Xk, logcr = dk$y)), 3))

cat("\n=== creatinine, same comparison, for reference ===\n")
cat("\ncovariates + log cystatin C:\n")
print(round(fit_families(dk$creat, cbind(Xk, logcys = log(dk$q))), 3))

## --- shape of each margin on the retained sample ------------------------
cat("\n--- shape diagnostics, eGFR >= 15 sample ---\n")
print(round(t(sapply(list(`cystatin C` = dk$q,
                          `log cystatin C` = log(dk$q),
                          `creatinine` = dk$creat,
                          `log creatinine` = log(dk$creat)),
                     function(x) c(skewness = skew(x), exc.kurt = kurt(x)))), 3))

cat("\nresiduals of log creatinine on covariates + log cystatin C:\n")
r <- resid(lm(dk$y ~ Xk + log(dk$q)))
cat(sprintf("  skew %.3f  exc.kurt %.3f  shapiro p %.3g\n",
            skew(r), kurt(r), shapiro.test(sample(r, 5000, replace = TRUE))$p.value))

## --- dependence vs kidney function, model-free ---------------------------
cat("\n--- Kendall's tau by eGFR stratum (stratifying on a covariate-based",
    "\n    staging variable, not on the two margins jointly) ---\n")
g <- cut(dk$egfr, c(15, 45, 60, 75, 90, 105, Inf),
         labels = c("15-45", "45-60", "60-75", "75-90", "90-105", ">105"))
print(round(sapply(levels(g), function(l) {
  ii <- which(g == l)
  c(n = length(ii), tau = cor(dk$q[ii], dk$y[ii], method = "kendall"))
}), 3))

cat("\n--- Kendall's tau by age quintile ---\n")
ga <- cut(dk$age, quantile(dk$age, seq(0, 1, length.out = 6)),
          include.lowest = TRUE, labels = paste0("A", 1:5))
print(round(sapply(levels(ga), function(l) {
  ii <- which(ga == l)
  c(n = length(ii), med.age = median(dk$age[ii]), tau = cor(dk$q[ii], dk$y[ii], method = "kendall"))
}), 3))

saveRDS(keep, "output/keep_final.rds")
