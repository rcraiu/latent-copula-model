## Build the NHANES 1999-2002 kidney-function analysis dataset.
##
## The latent variable W is (unobservable) renal impairment, the reciprocal
## scale of the glomerular filtration rate, for which the two continuous
## surrogates are two low-molecular-weight serum proteins:
##   Q = 1 / cystatin C            -- gamma margin, log link
##   Y = 1 / beta-2 microglobulin  -- normal margin
## Covariates: age, sex, BMI, systolic blood pressure.
##
## Both surrogates are on the reciprocal scale, and W is accordingly oriented as
## the latent filtration rate itself rather than as impairment: both markers
## increase with it. The reason is physiological before it is statistical. The
## serum concentration of a marker cleared by the kidneys is inversely
## proportional to the filtration rate, so it is the reciprocal of the
## concentration, not the concentration or its logarithm, that is linear in
## filtration; slopes of reciprocal serum creatinine against time are used
## clinically for precisely this reason, to track the progression of chronic
## kidney disease.
##
## The statistical case agrees, and it matters because the normal margin is the
## rigid one: Y given W is normal and W is normal, so the model implies a
## symmetric marginal distribution for Y and cannot do otherwise. Writing both
## margins as densities for the observed pair, so that specifications on
## different scales are comparable, the total AIC is 609 with Y on the log
## scale, -145 with Y reciprocal and Q left as a concentration, and -526 with
## both reciprocal. The probability integral transforms of the two margins have
## Kolmogorov-Smirnov distances from uniformity of 0.029 and 0.021 under the
## specification used here, against 0.055 and 0.054 on the log scale.
##
## This was not a cosmetic choice. Posterior predictive checks with Y on the log
## scale showed the model failing to reproduce the skewness of either margin,
## and the misfit biased the estimated dependence rather than merely the
## margins: it reversed the sign of the fitted trend in the copula parameter
## under the Gumbel family. See ppcheck.R.

## Both markers are freely filtered by the glomerulus and catabolised by the
## proximal tubule, so serum concentrations rise as filtration falls, and
## neither depends on muscle mass. That last point matters here: serum
## creatinine, the marker in routine clinical use, is generated in proportion to
## muscle mass and therefore depends on age, sex and body composition, which are
## the very quantities we wish to treat as covariates for W. Using two
## non-creatinine markers keeps the surrogates and the covariates distinct.
##
## Creatinine is deliberately *not* used as a surrogate. Besides the point
## above, it is reported on a 0.1 mg/dL grid in these cycles, giving only 31
## distinct values among several thousand observations with a sixth of the
## sample on a single value. A continuous margin cannot represent that, and with
## ties on that scale Sklar's theorem no longer identifies the copula. It is
## used here only to stage kidney disease for the exclusion below.

suppressMessages({
  library(haven)
  library(dplyr)
})

datadir <- "application-nhanes/data"
rx <- function(f) as.data.frame(read_xpt(file.path(datadir, paste0(f, ".xpt"))))

## --- the two surrogate markers, 1999-2004 surplus sera -------------------
## SSCARD_A carries the recalibrated cystatin C together with beta-2
## microglobulin. SSCYSTL and SSB2ML flag values outside the assay's reportable
## range, which are dropped rather than censored.
card <- rx("SSCARD_A") %>%
  select(SEQN, cystatinC = SSCYST, b2m = SSB2M,
         cyst_flag = SSCYSTL, b2m_flag = SSB2ML)

## --- covariates, 1999-2000 and 2001-2002 --------------------------------
demo <- bind_rows(
  rx("DEMO")   %>% select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1),
  rx("DEMO_B") %>% select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)
)

bmx <- bind_rows(
  rx("BMX")   %>% select(SEQN, BMXBMI),
  rx("BMX_B") %>% select(SEQN, BMXBMI)
)

bpx <- bind_rows(
  rx("BPX")   %>% select(SEQN, BPXSY1, BPXSY2, BPXSY3, BPXSY4),
  rx("BPX_B") %>% select(SEQN, BPXSY1, BPXSY2, BPXSY3, BPXSY4)
) %>%
  mutate(sbp = rowMeans(cbind(BPXSY1, BPXSY2, BPXSY3, BPXSY4), na.rm = TRUE)) %>%
  select(SEQN, sbp)

## serum creatinine, for staging only: LBXSCR in 1999-2000, LBDSCR in 2001-2002
cre <- bind_rows(
  rx("LAB18") %>% select(SEQN, creatinine = LBXSCR),
  rx("L40_B") %>% select(SEQN, creatinine = LBDSCR)
)

## --- merge and apply exclusions -----------------------------------------
## The inner join on demo restricts SSCARD_A to the 1999-2002 cycles.
raw <- card %>%
  inner_join(demo, by = "SEQN") %>%
  inner_join(bmx,  by = "SEQN") %>%
  inner_join(bpx,  by = "SEQN") %>%
  left_join(cre,   by = "SEQN")

cat("records in 1999-2002 with surplus sera measurements:", nrow(raw), "\n")

dat <- raw %>%
  filter(cyst_flag == 0, b2m_flag == 0) %>%
  filter(!is.na(cystatinC), !is.na(b2m), cystatinC > 0, b2m > 0) %>%
  filter(RIDAGEYR >= 20) %>%
  filter(!is.na(BMXBMI), !is.na(sbp), is.finite(sbp)) %>%
  transmute(
    SEQN,
    q     = 1 / cystatinC,            # gamma margin, reciprocal scale
    y     = 1 / b2m,                  # normal margin, reciprocal scale
    b2m       = b2m,
    cystatinC = cystatinC,
    creat = creatinine,
    age    = RIDAGEYR,
    female = as.numeric(RIAGENDR == 2),
    bmi    = BMXBMI,
    sbp    = sbp
  )

cat("adults with complete data:", nrow(dat), "\n")

## Exclude advanced chronic kidney disease, stages G3b to G5. Both markers lose
## their standing as clean filtration surrogates in this range: beta-2
## microglobulin in particular is raised by the inflammation and immune
## activation that accompany advanced disease, over and above any reduction in
## filtration, and participants in stage G5 are typically on dialysis, for which
## the relationship between serum markers and filtration breaks down altogether.
##
## Staging uses the creatinine-based CKD-EPI equation, so the exclusion is
## defined by a criterion external to both surrogates. This matters for the
## conditional model: creatinine is a noisy proxy for W itself rather than for
## the residual dependence between the surrogates, so restricting on it
## truncates the range of W over which the copula can be estimated but leaves
## the conditional laws of Q and Y given W, and the copula linking them,
## unaffected. Selecting instead on a function of the two surrogates would
## distort exactly the dependence the model sets out to estimate.
egfr_cr <- function(scr, age, female) {
  kap <- ifelse(female, 0.7, 0.9)
  alp <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / kap, 1)^alp * pmax(scr / kap, 1)^(-1.200) *
    0.9938^age * ifelse(female, 1.012, 1)
}
dat$egfr <- with(dat, egfr_cr(creat, age, female == 1))
n_adv <- sum(dat$egfr < 45, na.rm = TRUE)
dat <- dat[is.na(dat$egfr) | dat$egfr >= 45, ]
cat("excluded", n_adv, "in CKD stages G3b-G5 (eGFR < 45); final n =", nrow(dat), "\n")

## Standardized covariate matrix, no intercept (mu_q and mu_y are the intercepts)
X <- with(dat, cbind(age = age, female = female, bmi = bmi, sbp = sbp))
X <- scale(X)
attr(X, "scaled:center") -> xcenter
attr(X, "scaled:scale")  -> xscale
X <- matrix(as.numeric(X), nrow = nrow(dat), dimnames = list(NULL, colnames(X)))

NH <- list(q = dat$q, y = dat$y, X = X, dat = dat,
           xcenter = xcenter, xscale = xscale)

saveRDS(NH, "application-nhanes/nhanes_kidney.rds")

## --- descriptive summary ------------------------------------------------
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3

cat("\n--- marginal summaries ---\n")
print(round(rbind(
  `Q: 1/cystatin C`      = summary(dat$q),
  `B2M (mg/L)`           = summary(dat$b2m),
  `Y: 1/B2M`             = summary(dat$y)
), 3))

cat("\nskewness  Q = 1/cystatin C:", round(skew(dat$q), 3),
    " B2M:", round(skew(dat$b2m), 3),
    " Y = 1/B2M:", round(skew(dat$y), 3), "\n")

## The normal margin requires a symmetric residual; report the check.
ry <- residuals(lm(y ~ age + female + bmi + sbp, data = dat))
cat("residual of Y on covariates: skewness", round(skew(ry), 3),
    " Shapiro-Wilk", round(shapiro.test(sample(ry, min(4000, length(ry))))$statistic, 4), "\n")

## Granularity: the reason creatinine was replaced. Reported as the spacing
## between adjacent distinct values relative to the marginal standard
## deviation, and as the largest group of tied observations.
gran <- function(v, lab) {
  u <- sort(unique(v))
  cat(sprintf("%-12s distinct %5d   median spacing/sd %.4f   largest tie %5.2f%%\n",
              lab, length(u), median(diff(u)) / sd(v),
              100 * max(table(v)) / length(v)))
}
cat("\n--- granularity ---\n")
gran(dat$q, "Q: 1/cystatinC")
gran(dat$y, "Y: 1/B2M")
if (any(!is.na(dat$creat))) gran(log(dat$creat[!is.na(dat$creat)]), "log creat")

cat("\nPearson  cor(q, y) =", round(cor(dat$q, dat$y), 3),
    "  Kendall tau =", round(cor(dat$q, dat$y, method = "kendall"), 3),
    "  Spearman =", round(cor(dat$q, dat$y, method = "spearman"), 3), "\n")

ct <- cor.test(dat$q, dat$y, method = "kendall")
cat("Kendall test of independence: z =", round(ct$statistic, 2),
    " p =", format.pval(ct$p.value, digits = 3), "\n")

cat("\n--- covariates ---\n")
print(round(sapply(dat[, c("age", "female", "bmi", "sbp")], mean), 2))
print(round(sapply(dat[, c("age", "female", "bmi", "sbp")], sd), 2))
