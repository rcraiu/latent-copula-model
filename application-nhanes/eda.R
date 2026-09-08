## Exploratory analysis and figures for the NHANES kidney-function example.
##
##   Q = cystatin C (mg/L)             -- gamma margin, log link
##   Y = log beta-2 microglobulin      -- normal margin

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

base <- theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())

## ------------------------------------------------------------------------
## 1. Measurement granularity: why creatinine is not used as a surrogate
## ------------------------------------------------------------------------
## The model gives both surrogates continuous margins, and the copula is
## identified by Sklar's theorem only when the margins are continuous. Serum
## creatinine is reported to the nearest 0.1 mg/dL in these cycles, which is
## coarse relative to its own dispersion; the two markers used here are not.
gran <- function(v) {
  u <- sort(unique(v))
  c(distinct = length(u),
    spacing.sd = median(diff(u)) / sd(v),
    max.tie.pct = 100 * max(table(v)) / length(v))
}
cat("--- measurement granularity ---\n")
gtab <- rbind(
  `Q: cystatin C`           = gran(dat$q),
  `Y: log beta-2 microglob` = gran(dat$y),
  `log creatinine (unused)` = gran(log(dat$creat[!is.na(dat$creat)]))
)
print(round(gtab, 4))

## A figure making the contrast visible: the empirical CDF of each marker on
## its own standardised scale. A coarse grid shows as a staircase.
ecdf_df <- bind_rows(
  data.frame(z = scale(dat$q)[, 1],                     marker = "cystatin C"),
  data.frame(z = scale(dat$y)[, 1],                     marker = "log beta-2 microglobulin"),
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
## 2. Which transformation for which margin
## ------------------------------------------------------------------------
cand <- list(
  `cystatin C`               = dat$q,
  `log cystatin C`           = log(dat$q),
  `beta-2 microglobulin`     = dat$b2m,
  `log beta-2 microglobulin` = dat$y
)
tab <- t(sapply(cand, function(x)
  c(skewness = skew(x), excess.kurtosis = kurt(x), shapiro.W = sw(x))))
cat("\n--- marginal shape diagnostics ---\n"); print(round(tab, 3))

## The choice of margin is settled on the residual, not the raw marker, since
## the covariates account for part of the spread of each.
cat("\n--- residual shape after regression on the covariates ---\n")
rtab <- t(sapply(list(`log cystatin C`           = log(dat$q),
                      `log beta-2 microglobulin` = dat$y),
                 function(v) { r <- resid(lm(v ~ X))
                   c(skewness = skew(r), excess.kurtosis = kurt(r), shapiro.W = sw(r)) }))
print(round(rtab, 3))

## Gamma against lognormal for the margin taken to be gamma
g_gam <- glm(dat$q ~ X, family = Gamma(link = "log"))
cat("\ngamma GLM for cystatin C: dispersion =",
    round(summary(g_gam)$dispersion, 4),
    " => implied shape alpha ~", round(1 / summary(g_gam)$dispersion, 2), "\n")
## The lm on log(q) is a density for log(q), so its AIC needs the Jacobian
## 2*sum(log q) added before it can be set beside the gamma GLM.
cat("gamma GLM AIC =", round(AIC(g_gam), 1),
    "   lognormal AIC (same scale) =",
    round(AIC(lm(log(dat$q) ~ X)) + 2 * sum(log(dat$q)), 1), "\n")

## --- which marker takes which margin ------------------------------------
## Decided before fitting, by comparing the two assignments as densities for
## the pair (q, b2m) so that the totals are on a common scale, and by the
## calibration of each margin's probability integral transform.
assign_cmp <- function(gam, nor) {
  fg <- glm(gam ~ X, family = Gamma(link = "log"))
  fl <- lm(log(nor) ~ X)
  s  <- 1 / summary(fg)$dispersion
  ks <- function(u) unname(suppressWarnings(ks.test(u, "punif")$statistic))
  c(total.AIC = AIC(fg) + AIC(fl) + 2 * sum(log(nor)),
    KS.gamma  = ks(pgamma(gam, shape = s, scale = fitted(fg) / s)),
    KS.normal = ks(pnorm(resid(fl), sd = summary(fl)$sigma)))
}
cat("\n--- choice of margin (lower is better) ---\n")
print(round(rbind(
  `Q = cystatin C, Y = log B2M (used)` = assign_cmp(dat$q, dat$b2m),
  `Q = B2M, Y = log cystatin C`        = assign_cmp(dat$b2m, dat$q)
), 4))

## ------------------------------------------------------------------------
## 3. Joint behaviour
## ------------------------------------------------------------------------
u1 <- rank(dat$q) / (n + 1)
u2 <- rank(dat$y) / (n + 1)

cat("\n--- association ---\n")
cat("Pearson", round(cor(dat$q, dat$y), 3),
    " Kendall", round(cor(dat$q, dat$y, method = "kendall"), 3),
    " Spearman", round(cor(dat$q, dat$y, method = "spearman"), 3), "\n")

## Tail dependence, estimated nonparametrically, to see which families are
## plausible a priori. lambda_U(t) = P(U2 > t | U1 > t) and its lower analogue.
tailfun <- function(t) c(
  lower = mean(u2 <= t & u1 <= t) / t,
  upper = mean(u2 >  1 - t & u1 > 1 - t) / t
)
tt <- c(0.02, 0.05, 0.10, 0.15, 0.20)
cat("\n--- nonparametric tail dependence ---\n")
print(round(rbind(t = tt, sapply(tt, tailfun)), 3))

p1 <- ggplot(dat, aes(x = q, y = y)) +
  geom_point(alpha = 0.25, size = 0.7) +
  scale_x_log10() +
  labs(x = "cystatin C (mg/L, log scale)",
       y = "log beta-2 microglobulin (mg/L)") + base
ggsave("application-nhanes/figs/NHANES_scatter_orig.png", p1, width = 5.2, height = 4.4, dpi = 300)

p2 <- ggplot(data.frame(u1, u2), aes(x = u1, y = u2)) +
  geom_point(alpha = 0.25, size = 0.7) +
  labs(x = expression(hat(u)[1]~"(cystatin C)"),
       y = expression(hat(u)[2]~"(log beta-2 microglobulin)")) + base
ggsave("application-nhanes/figs/NHANES_scatter_pobs.png", p2, width = 5.2, height = 4.4, dpi = 300)

p3 <- ggplot(dat, aes(x = q)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "grey30", linewidth = 0.2) +
  coord_cartesian(xlim = c(0, 2)) +
  labs(x = "cystatin C (mg/L)", y = "count") + base
ggsave("application-nhanes/figs/NHANES_cystatinC_hist.png", p3, width = 5.2, height = 4.0, dpi = 300)

p4 <- ggplot(dat, aes(x = y)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "grey30", linewidth = 0.2) +
  labs(x = "log beta-2 microglobulin (mg/L)", y = "count") + base
ggsave("application-nhanes/figs/NHANES_logb2m_hist.png", p4, width = 5.2, height = 4.0, dpi = 300)

## ------------------------------------------------------------------------
## 4. Model-free look at whether the association varies with kidney function
## ------------------------------------------------------------------------
## We stratify on AGE, which is a covariate rather than a function of the two
## margins. Stratifying on the markers themselves induces a selection effect
## that attenuates the within-stratum association and is not interpretable.
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
