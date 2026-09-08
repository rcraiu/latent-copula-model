## Choose the gamma margin and the normal margin on the basis of fit,
## comparing candidate marginal families on a COMMON scale (density of the
## raw marker, so the log-Jacobian is included for log-scale fits).

suppressMessages({ library(dplyr) })

NH  <- readRDS("application-nhanes/nhanes_kidney.rds")
dat <- NH$dat
X   <- NH$X

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3

## how heavy is the tail?
cat("--- tail of the two markers ---\n")
print(round(quantile(dat$q,     c(.5, .9, .95, .99, .995, .999, 1)), 3))
print(round(quantile(dat$creat, c(.5, .9, .95, .99, .995, .999, 1)), 3))
cat("cystatin C > 3 mg/L:", sum(dat$q > 3), "   creatinine > 3 mg/dL:", sum(dat$creat > 3), "\n")

## ---- log-likelihood of the raw marker under each family, given covariates
loglik_gamma <- function(v, X) {
  fit <- glm(v ~ X, family = Gamma(link = "log"))
  mu  <- fitted(fit)
  nll <- function(la) -sum(dgamma(v, shape = exp(la), rate = exp(la) / mu, log = TRUE))
  la  <- optimize(nll, c(-3, 6))$minimum
  k   <- ncol(X) + 2
  list(loglik = -nll(la), alpha = exp(la), k = k, AIC = 2 * nll(la) + 2 * k)
}
loglik_lognormal <- function(v, X) {
  fit <- lm(log(v) ~ X)
  s   <- sqrt(mean(resid(fit)^2))
  ll  <- sum(dnorm(log(v), mean = fitted(fit), sd = s, log = TRUE)) - sum(log(v))
  k   <- ncol(X) + 2
  list(loglik = ll, k = k, AIC = -2 * ll + 2 * k)
}
loglik_normal <- function(v, X) {
  fit <- lm(v ~ X)
  s   <- sqrt(mean(resid(fit)^2))
  ll  <- sum(dnorm(v, mean = fitted(fit), sd = s, log = TRUE))
  k   <- ncol(X) + 2
  list(loglik = ll, k = k, AIC = -2 * ll + 2 * k)
}

compare <- function(v, label, X) {
  a <- loglik_gamma(v, X); b <- loglik_lognormal(v, X); d <- loglik_normal(v, X)
  cat(sprintf("\n%s:  gamma AIC %.1f (alpha=%.2f) | lognormal AIC %.1f | normal AIC %.1f\n",
              label, a$AIC, a$alpha, b$AIC, d$AIC))
  invisible(list(gamma = a, lognormal = b, normal = d))
}

cat("\n=== FULL SAMPLE (n =", nrow(dat), ") ===\n")
compare(dat$q,     "cystatin C", X)
compare(dat$creat, "creatinine", X)

## ---- restrict to the non-advanced-CKD range, the usual analytic sample --
keep <- dat$q <= 3 & dat$creat <= 3
cat("\n=== EXCLUDING cystatin C > 3 or creatinine > 3 (n =", sum(keep), ",",
    sum(!keep), "dropped) ===\n")
Xk <- X[keep, ]
compare(dat$q[keep],     "cystatin C", Xk)
compare(dat$creat[keep], "creatinine", Xk)

cat("\nshape diagnostics on restricted sample:\n")
print(round(t(sapply(list(`cystatin C` = dat$q[keep],
                          `log cystatin C` = log(dat$q[keep]),
                          `creatinine` = dat$creat[keep],
                          `log creatinine` = log(dat$creat[keep])),
                     function(x) c(skewness = skew(x), exc.kurt = kurt(x)))), 3))

cat("\nresidual shape, normal margin candidates (restricted sample):\n")
for (nm in c("log creatinine", "log cystatin C")) {
  v <- if (nm == "log creatinine") log(dat$creat[keep]) else log(dat$q[keep])
  r <- resid(lm(v ~ Xk))
  cat(sprintf("  %-16s skew %6.3f  exc.kurt %6.3f\n", nm, skew(r), kurt(r)))
}

## ---- model-free look at whether dependence varies with kidney function.
## Stratify on AGE (a covariate, not a function of the two margins) to avoid
## the selection bias induced by conditioning on the markers themselves.
cat("\n--- Kendall's tau by age quintile (restricted sample) ---\n")
dk  <- dat[keep, ]
grp <- cut(dk$age, quantile(dk$age, seq(0, 1, length.out = 6)),
           include.lowest = TRUE, labels = paste0("A", 1:5))
print(round(sapply(levels(grp), function(g) {
  ii <- which(grp == g); c(n = length(ii), tau = cor(dk$q[ii], dk$y[ii], method = "kendall"))
}), 3))

cat("\n--- Kendall's tau by strata of covariate-predicted marker level ---\n")
pred <- fitted(lm(scale(log(dk$q)) + scale(dk$y) ~ Xk))
grp2 <- cut(pred, quantile(pred, seq(0, 1, length.out = 6)),
            include.lowest = TRUE, labels = paste0("P", 1:5))
print(round(sapply(levels(grp2), function(g) {
  ii <- which(grp2 == g); c(n = length(ii), tau = cor(dk$q[ii], dk$y[ii], method = "kendall"))
}), 3))

saveRDS(keep, "output/keep_restricted.rds")
