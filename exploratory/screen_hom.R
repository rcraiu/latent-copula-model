## The insulin-resistance pairing, with the clinically standard exclusion.
##
## Latent W = insulin resistance. It cannot be measured directly: the gold
## standard is the euglycaemic-hyperinsulinaemic clamp, infeasible in a survey.
## Fasting insulin and fasting glucose are the two standard surrogates, which
## practice combines into the ad hoc product HOMA-IR = glucose x insulin / 405.
##
## Diagnosed diabetics and those on glucose-lowering medication are excluded,
## as is standard when studying insulin resistance in the general population,
## because treatment breaks the link between the surrogates and the latent
## quantity. This is a clinical exclusion, not one chosen to improve fit.

suppressMessages({ library(haven); library(dplyr) })
rx <- function(f) as.data.frame(read_xpt(file.path("data", paste0(f, ".xpt"))))

demo <- rx("DEMO_J") %>% select(SEQN, age = RIDAGEYR, sex = RIAGENDR,
                                eth = RIDRETH3)
ins  <- rx("INS_J")  %>% select(SEQN, insulin = LBXIN)
glu  <- rx("GLU_J")  %>% select(SEQN, glucose = LBXGLU)
bmx  <- rx("BMX_J")  %>% select(SEQN, bmi = BMXBMI, waist = BMXWAIST)
diq  <- rx("DIQ_J")  %>% select(SEQN, dx = DIQ010, insulin_rx = DIQ050,
                                pills = DIQ070)
ghb  <- rx("GHB_J")  %>% select(SEQN, a1c = LBXGH)

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3

D <- demo %>% filter(age >= 20) %>%
  left_join(ins, "SEQN") %>% left_join(glu, "SEQN") %>%
  left_join(bmx, "SEQN") %>% left_join(diq, "SEQN") %>% left_join(ghb, "SEQN") %>%
  mutate(female = as.numeric(sex == 2))

step <- function(d, lab) { cat(sprintf("%-46s n = %d\n", lab, nrow(d))); d }

d <- D %>% filter(!is.na(insulin), !is.na(glucose), insulin > 0, glucose > 0) %>%
  step("with fasting insulin and glucose") %>%
  filter(!is.na(bmi), !is.na(waist)) %>% step("with body measures") %>%
  filter(!is.na(dx)) %>% step("with diabetes questionnaire")

d_all <- d
d <- d %>% filter(dx != 1) %>% step("excluding diagnosed diabetes") %>%
  filter(is.na(insulin_rx) | insulin_rx != 1) %>%
  filter(is.na(pills) | pills != 1) %>% step("excluding glucose-lowering medication") %>%
  filter(!is.na(a1c), a1c < 6.5) %>% step("excluding undiagnosed diabetes (HbA1c >= 6.5)")

gran <- function(v, nm) {
  u <- sort(unique(v))
  cat(sprintf("%-22s distinct = %5d   max tie = %5.2f%%   min gap = %.3g\n",
              nm, length(u), 100 * max(table(v)) / length(v), min(diff(u))))
}
cat("\n--- granularity ---\n")
gran(d$insulin, "fasting insulin (uU/mL)")
gran(d$glucose, "fasting glucose (mg/dL)")

cat("\n--- marginal shape ---\n")
print(round(t(sapply(list(insulin = d$insulin, `log insulin` = log(d$insulin),
                          glucose = d$glucose, `log glucose` = log(d$glucose)),
                     function(x) c(skew = skew(x), exc.kurt = kurt(x)))), 3))

Z <- scale(as.matrix(d[, c("age", "female", "bmi", "waist")]))
Z <- matrix(as.numeric(Z), nrow(d))

## gamma vs lognormal for insulin, on the raw-marker scale
gl  <- glm(d$insulin ~ Z, family = Gamma(link = "log")); mu <- fitted(gl)
nll <- function(la) -sum(dgamma(d$insulin, shape = exp(la), rate = exp(la) / mu, log = TRUE))
la  <- optimize(nll, c(-4, 8))$minimum
k   <- ncol(Z) + 2
gA  <- 2 * nll(la) + 2 * k
l1  <- lm(log(d$insulin) ~ Z); s <- sqrt(mean(resid(l1)^2))
lA  <- -2 * (sum(dnorm(log(d$insulin), fitted(l1), s, log = TRUE)) - sum(log(d$insulin))) + 2 * k
cat(sprintf("\ninsulin: gamma AIC %.1f (alpha %.2f) | lognormal AIC %.1f | per-obs %+.4f\n",
            gA, exp(la), lA, (gA - lA) / nrow(d)))

## same comparison conditioning also on the other surrogate, a proxy for W
Z2 <- cbind(Z, logglu = log(d$glucose))
gl2 <- glm(d$insulin ~ Z2, family = Gamma(link = "log")); mu2 <- fitted(gl2)
nll2 <- function(la) -sum(dgamma(d$insulin, shape = exp(la), rate = exp(la) / mu2, log = TRUE))
la2 <- optimize(nll2, c(-4, 8))$minimum
k2  <- ncol(Z2) + 2
gA2 <- 2 * nll2(la2) + 2 * k2
l2  <- lm(log(d$insulin) ~ Z2); s2 <- sqrt(mean(resid(l2)^2))
lA2 <- -2 * (sum(dnorm(log(d$insulin), fitted(l2), s2, log = TRUE)) - sum(log(d$insulin))) + 2 * k2
cat(sprintf("insulin | covariates + log glucose: gamma AIC %.1f (alpha %.2f) | lognormal AIC %.1f | per-obs %+.4f\n",
            gA2, exp(la2), lA2, (gA2 - lA2) / nrow(d)))

cat("\n--- normality of log glucose ---\n")
r <- resid(lm(log(d$glucose) ~ Z))
cat(sprintf("  residual skew %+.3f  exc.kurt %+.3f\n", skew(r), kurt(r)))
r2 <- resid(lm(log(d$glucose) ~ Z2[, 1:4] + log(d$insulin)))
cat(sprintf("  residual skew %+.3f  exc.kurt %+.3f  (also given log insulin)\n",
            skew(r2), kurt(r2)))

cat("\n--- dependence ---\n")
ct <- cor.test(d$insulin, d$glucose, method = "kendall")
cat("Kendall tau =", round(ct$estimate, 3), " z =", round(ct$statistic, 2),
    " p =", format.pval(ct$p.value, digits = 3), "\n")
cat("Spearman =", round(cor(d$insulin, d$glucose, method = "spearman"), 3), "\n")

cat("\n--- Kendall's tau by BMI quintile (beta-cell compensation) ---\n")
for (v in c("bmi", "age", "waist")) {
  g <- cut(d[[v]], quantile(d[[v]], seq(0, 1, length.out = 6)),
           include.lowest = TRUE, labels = paste0("Q", 1:5))
  taus <- sapply(levels(g), function(l) {
    ii <- which(g == l); cor(d$insulin[ii], d$glucose[ii], method = "kendall")
  })
  cat(sprintf("  by %-6s %s  (spread %.3f)\n", v,
              paste(sprintf("%+.3f", taus), collapse = "  "), diff(range(taus))))
}

saveRDS(d, "application-nhanes/nhanes_homa.rds")
