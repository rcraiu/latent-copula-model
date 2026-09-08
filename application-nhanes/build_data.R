## Build the NHANES 1999-2002 kidney-function analysis dataset.
##
## The latent variable W is the (unobservable) true glomerular filtration rate,
## for which the two continuous surrogates are
##   Q = serum cystatin C (mg/L)        -- gamma margin
##   Y = log serum creatinine (mg/dL)   -- normal margin
## Covariates: age, sex, BMI, systolic blood pressure.

suppressMessages({
  library(haven)
  library(dplyr)
})

datadir <- "application-nhanes/data"
rx <- function(f) as.data.frame(read_xpt(file.path(datadir, paste0(f, ".xpt"))))

## --- cystatin C (surplus sera), 1999-2000 and 2001-2002 -----------------
cys <- bind_rows(
  rx("SSCYST_A") %>% select(SEQN, SSCYPC),
  rx("SSCYST_B") %>% select(SEQN, SSCYPC)
) %>% rename(cystatinC = SSCYPC)

## --- serum creatinine: LBXSCR in 1999-2000, LBDSCR in 2001-2002 ---------
cre <- bind_rows(
  rx("LAB18") %>% select(SEQN, creatinine = LBXSCR),
  rx("L40_B") %>% select(SEQN, creatinine = LBDSCR)
)

## --- covariates ---------------------------------------------------------
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

## --- merge and apply exclusions -----------------------------------------
raw <- cys %>%
  inner_join(cre, by = "SEQN") %>%
  inner_join(demo, by = "SEQN") %>%
  inner_join(bmx, by = "SEQN") %>%
  inner_join(bpx, by = "SEQN")

cat("merged records with both markers:", nrow(raw), "\n")

dat <- raw %>%
  filter(RIDAGEYR >= 20) %>%
  filter(!is.na(cystatinC), !is.na(creatinine), cystatinC > 0, creatinine > 0) %>%
  filter(!is.na(BMXBMI), !is.na(sbp), is.finite(sbp)) %>%
  transmute(
    SEQN,
    q      = cystatinC,               # gamma margin
    y      = log(creatinine),         # normal margin
    creat  = creatinine,
    age    = RIDAGEYR,
    female = as.numeric(RIAGENDR == 2),
    bmi    = BMXBMI,
    sbp    = sbp
  )

cat("analysis sample after exclusions:", nrow(dat), "\n")

## Exclude kidney failure (CKD stage G5). This is a clinical exclusion, not a
## marker cutoff chosen for fit: such participants are typically on dialysis,
## for which the relationship between filtration markers and GFR breaks down.
## eGFR is used only for staging here; the surrogates are the raw markers, so
## that age and sex enter the model as covariates for W rather than being
## baked into the response as they are in the CKD-EPI equations.
egfr_cr <- function(scr, age, female) {
  kap <- ifelse(female, 0.7, 0.9)
  alp <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / kap, 1)^alp * pmax(scr / kap, 1)^(-1.200) *
    0.9938^age * ifelse(female, 1.012, 1)
}
dat$egfr <- egfr_cr(dat$creat, dat$age, dat$female == 1)
n_g5 <- sum(dat$egfr < 15)
dat <- dat[dat$egfr >= 15, ]
cat("excluded", n_g5, "with eGFR < 15; final n =", nrow(dat), "\n")

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
cat("\n--- marginal summaries ---\n")
print(round(rbind(
  `cystatin C (mg/L)` = summary(dat$q),
  `creatinine (mg/dL)` = summary(dat$creat),
  `log creatinine` = summary(dat$y)
), 3))

cat("\nskewness  cystatin C:", round(mean((dat$q - mean(dat$q))^3) / sd(dat$q)^3, 3),
    " creatinine:", round(mean((dat$creat - mean(dat$creat))^3) / sd(dat$creat)^3, 3),
    " log creatinine:", round(mean((dat$y - mean(dat$y))^3) / sd(dat$y)^3, 3), "\n")

cat("\nPearson  cor(q, y) =", round(cor(dat$q, dat$y), 3),
    "  Kendall tau =", round(cor(dat$q, dat$y, method = "kendall"), 3),
    "  Spearman =", round(cor(dat$q, dat$y, method = "spearman"), 3), "\n")

ct <- cor.test(dat$q, dat$y, method = "kendall")
cat("Kendall test of independence: z =", round(ct$statistic, 2),
    " p =", format.pval(ct$p.value, digits = 3), "\n")

cat("\n--- covariates ---\n")
print(round(sapply(dat[, c("age", "female", "bmi", "sbp")], mean), 2))
print(round(sapply(dat[, c("age", "female", "bmi", "sbp")], sd), 2))
