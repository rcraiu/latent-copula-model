## NHANES serum creatinine is reported on a 0.1 mg/dL grid in every cycle
## (LLOD 0.1), giving only 31 distinct values here, so it cannot serve as a
## continuous margin. Candidate replacement partner for cystatin C: the urine
## albumin-to-creatinine ratio (ACR), the other axis of KDIGO CKD staging.

suppressMessages({ library(haven); library(dplyr) })
rx <- function(f) as.data.frame(read_xpt(file.path("data", paste0(f, ".xpt"))))

alb <- bind_rows(rx("LAB16") %>% select(SEQN, URXUMA, URXUCR),
                 rx("L16_B") %>% select(SEQN, URXUMA, URXUCR))
NH  <- readRDS("application-nhanes/nhanes_kidney.rds")

d <- NH$dat %>% inner_join(alb, by = "SEQN") %>%
  filter(!is.na(URXUMA), !is.na(URXUCR), URXUMA > 0, URXUCR > 0) %>%
  mutate(acr = URXUMA / URXUCR * 100)      # mg albumin per g creatinine

cat("n with cystatin C, ACR and covariates:", nrow(d), "\n\n")

gran <- function(v, nm) {
  u <- sort(unique(v))
  cat(sprintf("%-22s distinct=%5d  maxtie=%5.2f%%  min gap=%.4g\n",
              nm, length(u), 100 * max(table(v)) / length(v), min(diff(u))))
}
gran(d$q,      "cystatin C (mg/L)")
gran(d$URXUMA, "urine albumin (mg/L)")
gran(d$URXUCR, "urine creatinine (mg/dL)")
gran(d$acr,    "ACR (mg/g)")
gran(d$creat,  "serum creatinine (mg/dL)")

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3
cat("\n--- marginal shape ---\n")
print(round(t(sapply(list(`cystatin C` = d$q, `log cystatin C` = log(d$q),
                          `ACR` = d$acr, `log ACR` = log(d$acr),
                          `urine albumin` = d$URXUMA,
                          `log urine albumin` = log(d$URXUMA)),
                     function(x) c(skew = skew(x), exc.kurt = kurt(x)))), 3))

## proper family comparison on a common scale, conditioning on covariates
X <- NH$X[match(d$SEQN, NH$dat$SEQN), ]
fam_cmp <- function(v, Z, nm) {
  gl <- glm(v ~ Z, family = Gamma(link = "log")); mu <- fitted(gl)
  nll <- function(la) -sum(dgamma(v, shape = exp(la), rate = exp(la) / mu, log = TRUE))
  la <- optimize(nll, c(-4, 8))$minimum
  k  <- ncol(as.matrix(Z)) + 2
  gA <- 2 * nll(la) + 2 * k
  l1 <- lm(log(v) ~ Z); s <- sqrt(mean(resid(l1)^2))
  ll <- sum(dnorm(log(v), fitted(l1), s, log = TRUE)) - sum(log(v))
  cat(sprintf("%-16s gamma AIC %9.1f (alpha %6.2f) | lognormal AIC %9.1f | per-obs %+.4f\n",
              nm, gA, exp(la), -2 * ll + 2 * k, (gA - (-2 * ll + 2 * k)) / length(v)))
}
cat("\n--- marginal family comparison, given covariates ---\n")
fam_cmp(d$q,   X, "cystatin C")
fam_cmp(d$acr, X, "ACR")

cat("\n--- residual normality of candidate NORMAL margins ---\n")
for (nm in c("log cystatin C", "log ACR")) {
  v <- if (nm == "log cystatin C") log(d$q) else log(d$acr)
  r <- resid(lm(v ~ X))
  cat(sprintf("  %-16s skew %+.3f  exc.kurt %+.3f\n", nm, skew(r), kurt(r)))
}

cat("\n--- dependence between cystatin C and ACR ---\n")
cat("Kendall tau =", round(cor(d$q, d$acr, method = "kendall"), 3),
    " Spearman =", round(cor(d$q, d$acr, method = "spearman"), 3), "\n")
ct <- cor.test(d$q, d$acr, method = "kendall")
cat("test of independence: z =", round(ct$statistic, 2),
    " p =", format.pval(ct$p.value, digits = 3), "\n")

cat("\n--- Kendall's tau by age quintile (cystatin C vs ACR) ---\n")
g <- cut(d$age, quantile(d$age, seq(0, 1, length.out = 6)),
         include.lowest = TRUE, labels = paste0("A", 1:5))
print(round(sapply(levels(g), function(l) {
  ii <- which(g == l)
  c(n = length(ii), med.age = median(d$age[ii]),
    tau = cor(d$q[ii], d$acr[ii], method = "kendall"))
}), 3))

saveRDS(d, "nhanes_acr_candidate.rds")
