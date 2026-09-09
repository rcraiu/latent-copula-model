## Exploratory analysis and figures for the NHANES kidney-function example.
##
##   Q = 1 / cystatin C            -- gamma margin, log link
##   Y = 1 / beta-2 microglobulin  -- normal margin
##
## Both surrogates are on the reciprocal, or clearance, scale, on which the
## concentration of a renally cleared marker is linear in the filtration rate.
## W is therefore oriented as latent filtration: both markers increase with it.

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

NH  <- readRDS("application-nhanes/nhanes_kidney.rds")
dat <- NH$dat
X   <- NH$X
n   <- nrow(dat)
dir.create("application-nhanes/figs", recursive = TRUE, showWarnings = FALSE)

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3
sw   <- function(x) unname(shapiro.test(x[sample(length(x), min(4000, length(x)))])$statistic)
ks   <- function(u) unname(suppressWarnings(ks.test(u, "punif")$statistic))

base <- theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())

## ------------------------------------------------------------------------
## 1. Measurement granularity: why creatinine is not used as a surrogate
## ------------------------------------------------------------------------
## The model gives both surrogates continuous margins, and Sklar's theorem
## identifies the copula uniquely only when the margins are continuous. Serum
## creatinine is reported to the nearest 0.1 mg/dL in these cycles, which is
## coarse relative to its own dispersion; the two markers used here are not.
gran <- function(v) {
  u <- sort(unique(v))
  c(distinct = length(u), spacing.sd = median(diff(u)) / sd(v),
    max.tie.pct = 100 * max(table(v)) / length(v))
}
cat("--- measurement granularity ---\n")
print(round(rbind(
  `Q: 1/cystatin C`         = gran(dat$q),
  `Y: 1/B2M`                = gran(dat$y),
  `log creatinine (unused)` = gran(log(dat$creat[!is.na(dat$creat)]))
), 4))

ecdf_df <- bind_rows(
  data.frame(z = scale(dat$q)[, 1], marker = "1 / cystatin C"),
  data.frame(z = scale(dat$y)[, 1], marker = "1 / beta-2 microglobulin"),
  data.frame(z = scale(log(dat$creat[!is.na(dat$creat)]))[, 1],
             marker = "log creatinine (not used)")
)
p0 <- ggplot(ecdf_df, aes(x = z)) +
  stat_ecdf(geom = "step", linewidth = 0.5) +
  facet_wrap(~ marker, nrow = 1) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  labs(x = "standardised marker", y = "empirical distribution function") +
  base + theme(strip.background = element_rect(fill = "grey92"))
ggsave("application-nhanes/figs/NHANES_granularity.png", p0, width = 9.6, height = 3.4, dpi = 300)

## ------------------------------------------------------------------------
## 2. Choice of scale and of margin
## ------------------------------------------------------------------------
## The normal margin is the rigid one: Y given W is normal and W is normal, so
## the model implies a symmetric marginal Y and cannot do otherwise. The scale
## for Y therefore has to be one on which the residual is symmetric. The gamma
## margin has a shape parameter and is more forgiving, but is not indifferent.
cat("\n--- marginal shape on candidate scales ---\n")
cand <- list(
  `cystatin C`       = dat$cystatinC,  `log cystatin C` = log(dat$cystatinC),
  `1/cystatin C`     = dat$q,
  `B2M`              = dat$b2m,        `log B2M`        = log(dat$b2m),
  `1/B2M`            = dat$y
)
print(round(t(sapply(cand, function(x)
  c(skewness = skew(x), excess.kurtosis = kurt(x), shapiro.W = sw(x)))), 3))

cat("\n--- residual shape after regression on the covariates ---\n")
print(round(t(sapply(cand, function(v) { r <- resid(lm(v ~ X))
  c(skewness = skew(r), excess.kurtosis = kurt(r), shapiro.W = sw(r)) })), 3))

## The Box-Cox power that symmetrises each marker's residual, for comparison
## with the reciprocal. Both land close to -1, which is the clearance scale.
cat("\n--- Box-Cox power symmetrising the residual (reciprocal is -1) ---\n")
for (nm in c("cystatinC", "b2m")) {
  v <- dat[[nm]]
  f <- function(l) abs(skew(resid(lm(if (abs(l) < 1e-8) log(v) else (v^l - 1)/l ~ X))))
  cat(sprintf("  %-12s optimum %6.3f   |skew| there %.4f\n",
              nm, optimize(f, c(-2.5, 1))$minimum, optimize(f, c(-2.5, 1))$objective))
}

## Total AIC of the two marginal fits, all written as densities for the observed
## pair (cystatin C, B2M) so that specifications on different scales compare.
spec <- function(gam, jg, nor, jn, lab) {
  fg <- glm(gam ~ X, family = Gamma(link = "log")); s <- 1 / summary(fg)$dispersion
  fl <- lm(nor ~ X)
  cat(sprintf("%-42s %10.1f %9.4f %9.4f\n", lab,
      AIC(fg) - 2*jg + AIC(fl) - 2*jn,
      ks(pgamma(gam, shape = s, scale = fitted(fg)/s)),
      ks(pnorm(resid(fl), sd = summary(fl)$sigma))))
}
cy <- dat$cystatinC; b <- dat$b2m
cat("\n--- choice of specification (lower AIC and KS are better) ---\n")
cat(sprintf("%-42s %10s %9s %9s\n", "specification", "total AIC", "KS gamma", "KS normal"))
spec(cy, 0, log(b), sum(log(1/b)), "Q=cystatin C,   Y=log B2M")
spec(b, 0, log(cy), sum(log(1/cy)), "Q=B2M,          Y=log cystatin C")
spec(cy, 0, 1/b, sum(log(1/b^2)), "Q=cystatin C,   Y=1/B2M")
spec(1/cy, sum(log(1/cy^2)), 1/b, sum(log(1/b^2)), "Q=1/cystatin C, Y=1/B2M   (used)")
spec(1/b, sum(log(1/b^2)), 1/cy, sum(log(1/cy^2)), "Q=1/B2M,        Y=1/cystatin C")

## ------------------------------------------------------------------------
## 3. Joint behaviour
## ------------------------------------------------------------------------
u1 <- rank(dat$q) / (n + 1)
u2 <- rank(dat$y) / (n + 1)

cat("\n--- association ---\n")
cat("Pearson", round(cor(dat$q, dat$y), 3),
    " Kendall", round(cor(dat$q, dat$y, method = "kendall"), 3),
    " Spearman", round(cor(dat$q, dat$y, method = "spearman"), 3), "\n")

## Tail dependence, to see which families are plausible a priori. Note the
## orientation: on the clearance scale the *lower* tail is poor filtration.
tailfun <- function(t) c(lower = mean(u2 <= t & u1 <= t) / t,
                         upper = mean(u2 > 1 - t & u1 > 1 - t) / t)
tt <- c(0.02, 0.05, 0.10, 0.15, 0.20)
cat("\n--- nonparametric tail dependence (lower tail = poor filtration) ---\n")
print(round(rbind(t = tt, sapply(tt, tailfun)), 3))

p1 <- ggplot(dat, aes(x = q, y = y)) +
  geom_point(alpha = 0.25, size = 0.7) +
  labs(x = "1 / cystatin C  (L/mg)", y = "1 / beta-2 microglobulin  (L/mg)") + base
ggsave("application-nhanes/figs/NHANES_scatter_orig.png", p1, width = 5.2, height = 4.4, dpi = 300)

p2 <- ggplot(data.frame(u1, u2), aes(x = u1, y = u2)) +
  geom_point(alpha = 0.25, size = 0.7) +
  labs(x = expression(hat(u)[1]~"(1 / cystatin C)"),
       y = expression(hat(u)[2]~"(1 / beta-2 microglobulin)")) + base
ggsave("application-nhanes/figs/NHANES_scatter_pobs.png", p2, width = 5.2, height = 4.4, dpi = 300)

p3 <- ggplot(dat, aes(x = q)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "grey30", linewidth = 0.2) +
  labs(x = "1 / cystatin C  (L/mg)", y = "count") + base
ggsave("application-nhanes/figs/NHANES_cystatinC_hist.png", p3, width = 5.2, height = 4.0, dpi = 300)

p4 <- ggplot(dat, aes(x = y)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "grey30", linewidth = 0.2) +
  labs(x = "1 / beta-2 microglobulin  (L/mg)", y = "count") + base
ggsave("application-nhanes/figs/NHANES_logb2m_hist.png", p4, width = 5.2, height = 4.0, dpi = 300)

## ------------------------------------------------------------------------
## 4. Does the association vary with kidney function?
## ------------------------------------------------------------------------
## Stratify on age, a covariate, rather than on any function of the two margins:
## conditioning on the latter induces a selection effect which attenuates the
## within-stratum association and is not interpretable.
grp <- cut(dat$age, quantile(dat$age, seq(0, 1, length.out = 6)),
           include.lowest = TRUE, labels = paste0("A", 1:5))
loc <- sapply(levels(grp), function(g) {
  ii <- which(grp == g)
  ct <- cor.test(dat$q[ii], dat$y[ii], method = "kendall")
  c(n = length(ii), med.age = median(dat$age[ii]),
    tau = unname(ct$estimate), p = ct$p.value)
})
cat("\n--- Kendall's tau by age quintile ---\n")
print(round(loc, 4))

p5 <- ggplot(data.frame(age = as.numeric(loc["med.age", ]),
                        tau = as.numeric(loc["tau", ])),
             aes(x = age, y = tau)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.4) +
  labs(x = "median age within quintile (years)",
       y = expression("Kendall's "*tau*" between the two surrogates")) + base
ggsave("application-nhanes/figs/NHANES_tau_by_age.png", p5, width = 5.2, height = 4.0, dpi = 300)

cat("\nfigures written to application-nhanes/figs/\n")
