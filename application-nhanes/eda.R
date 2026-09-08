## Exploratory analysis and figures for the NHANES kidney-function example.

suppressMessages({
  library(dplyr)
  library(ggplot2)
})

NH  <- readRDS("application-nhanes/nhanes_kidney.rds")
dat <- NH$dat
dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3

## --- which transformation for which margin? -----------------------------
cand <- list(
  `cystatin C`      = dat$q,
  `log cystatin C`  = log(dat$q),
  `creatinine`      = dat$creat,
  `log creatinine`  = log(dat$creat)
)
tab <- t(sapply(cand, function(x) c(skewness = skew(x), excess.kurtosis = kurt(x),
                                    shapiro.W = shapiro.test(sample(x, 5000, replace = TRUE))$statistic)))
cat("--- marginal shape diagnostics ---\n"); print(round(tab, 3))

## Gamma vs lognormal fit to the *conditional* margin, via GLMs on covariates
X <- NH$X
g_gam <- glm(dat$q ~ X, family = Gamma(link = "log"))
cat("\ngamma GLM for cystatin C: dispersion =", round(summary(g_gam)$dispersion, 4),
    " => implied shape alpha ~", round(1 / summary(g_gam)$dispersion, 2), "\n")
cat("gamma GLM AIC =", round(AIC(g_gam), 1), "\n")
g_ln <- lm(log(dat$q) ~ X)
cat("lognormal (lm on log q) AIC =", round(AIC(g_ln), 1), "\n")

l_norm <- lm(dat$y ~ X)
cat("\nnormal lm for log creatinine: residual skewness =", round(skew(resid(l_norm)), 3),
    " excess kurtosis =", round(kurt(resid(l_norm)), 3), "\n")

## --- pseudo-observations ------------------------------------------------
n <- nrow(dat)
u1 <- rank(dat$q) / (n + 1)
u2 <- rank(dat$y) / (n + 1)

## --- figures ------------------------------------------------------------
base <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

p1 <- ggplot(dat, aes(x = q, y = y)) +
  geom_point(alpha = 0.25, size = 0.7) +
  scale_x_log10() +
  labs(x = "cystatin C (mg/L, log scale)", y = "log serum creatinine (mg/dL)") +
  base
ggsave("output/figs/NHANES_scatter_orig.png", p1, width = 5.2, height = 4.4, dpi = 300)

p2 <- ggplot(data.frame(u1, u2), aes(x = u1, y = u2)) +
  geom_point(alpha = 0.25, size = 0.7) +
  labs(x = expression(hat(u)[1]~"(cystatin C)"), y = expression(hat(u)[2]~"(log creatinine)")) +
  base
ggsave("output/figs/NHANES_scatter_pobs.png", p2, width = 5.2, height = 4.4, dpi = 300)

p3 <- ggplot(dat, aes(x = q)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "grey30", linewidth = 0.2) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(x = "cystatin C (mg/L)", y = "count") + base
ggsave("output/figs/NHANES_cystatinC_hist.png", p3, width = 5.2, height = 4.0, dpi = 300)

p4 <- ggplot(dat, aes(x = y)) +
  geom_histogram(bins = 60, fill = "grey70", colour = "grey30", linewidth = 0.2) +
  labs(x = "log serum creatinine (mg/dL)", y = "count") + base
ggsave("output/figs/NHANES_logcreat_hist.png", p4, width = 5.2, height = 4.0, dpi = 300)

## Model-free look at whether the association varies with kidney function.
## We stratify on AGE, which is a covariate rather than a function of the two
## margins. Stratifying on the markers themselves (or on eGFR, which is a
## deterministic function of creatinine) induces a selection effect that
## attenuates the within-stratum association and is not interpretable here.
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
ggsave("output/figs/NHANES_tau_by_age.png", p5, width = 5.2, height = 4.0, dpi = 300)

cat("\nfigures written to figs/\n")
