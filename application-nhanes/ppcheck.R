## Posterior predictive check for the NHANES kidney fit.
##
## The margin diagnostics in eda.R are approximations: they fit a gamma GLM and
## a linear model on the covariates alone, ignoring the latent variable, and so
## cannot say whether the *model's* margins are adequate. The margin of Q under
## the model is a gamma mixed over W, which is more skewed than any gamma, and
## the margin of Y is a normal mixed over W. This script therefore draws
## replicate data sets from the posterior predictive distribution and compares
## features of the observed data with their predictive distributions.
##
## Usage: Rscript ppcheck.R [family] [nrep]

suppressMessages(library(ggplot2))
source("R/copula_densities.R")
suppressMessages(library(VineCopula))

family <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "Frank"
R      <- if (length(commandArgs(TRUE)) >= 2) as.integer(commandArgs(TRUE)[2]) else 500

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X; n <- length(q)

VC <- c(Gaussian = 1, Clayton = 3, Gumbel = 4, Frank = 5, Joe = 6,
        SurvClayton = 13)

fit <- readRDS(sprintf("output/fits/%s_chain1.rds", family))
S   <- fit$PAR[, (fit$burn + 1):fit$M, drop = FALSE]
dr  <- sample(ncol(S), R)

## Features to check: shape of each margin, the strength of the association,
## and its asymmetry between the tails, which is what discriminates families.
feat <- function(qq, yy) {
  u1 <- rank(qq) / (length(qq) + 1); u2 <- rank(yy) / (length(yy) + 1)
  skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
  lam <- function(t, up) if (up) mean(u1 > 1 - t & u2 > 1 - t) / t else
                                 mean(u1 <= t & u2 <= t) / t
  c(`skew q` = skew(qq), `sd q` = sd(qq), `q97.5` = quantile(qq, 0.975, names = FALSE),
    `skew y` = skew(yy), `sd y` = sd(yy), `y97.5` = quantile(yy, 0.975, names = FALSE),
    `kendall` = cor(qq, yy, method = "kendall"),
    `lambda lo` = lam(0.05, FALSE), `lambda up` = lam(0.05, TRUE),
    `tail asym` = lam(0.05, TRUE) - lam(0.05, FALSE))
}

obs <- feat(q, y)
rep <- matrix(NA_real_, R, length(obs), dimnames = list(NULL, names(obs)))

for (r in seq_len(R)) {
  p <- S[, dr[r]]
  sW <- if ("s2W" %in% names(p)) sqrt(p[["s2W"]]) else 1
  W <- as.numeric(X %*% p[c("beta1", "beta2", "beta3", "beta4")]) + rnorm(n, sd = sW)
  al <- exp(p[["lalpha"]]); sg <- exp(p[["lsigma"]])
  mq <- exp(p[["muq"]] + p[["lamq"]] * W)
  my <- p[["muy"]] + p[["lamy"]] * W
  if (family == "Independence") {
    u1 <- runif(n); u2 <- runif(n)
  } else {
    th <- theta_of_w(W, c(p[["xi0"]], p[["xi1"]]), family, 0.1)
    u1 <- runif(n)
    ## conditional sampling: U2 | U1 = u1 from the h-inverse of the copula
    u2 <- suppressWarnings(BiCopHinv1(u1, runif(n), family = VC[[family]], par = th))
  }
  rep[r, ] <- feat(qgamma(u1, shape = al, scale = mq / al),
                   qnorm(u2, mean = my, sd = sg))
}

pp <- colMeans(sweep(rep, 2, obs, ">="))
out <- data.frame(feature = names(obs), observed = obs,
                  pred.mean = colMeans(rep),
                  pred.lo = apply(rep, 2, quantile, 0.025),
                  pred.hi = apply(rep, 2, quantile, 0.975),
                  post.p = pp)
cat(sprintf("posterior predictive check, %s copula, %d replicates\n\n", family, R))
print(out, row.names = FALSE, digits = 3)
cat("\npost.p is Pr(replicate >= observed); values near 0 or 1 indicate misfit.\n")

dir.create("application-nhanes/figs", recursive = TRUE, showWarnings = FALSE)
long <- do.call(rbind, lapply(names(obs), function(nm)
  data.frame(feature = nm, value = rep[, nm], observed = obs[[nm]])))
p <- ggplot(long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey75", colour = "grey35", linewidth = 0.15) +
  geom_vline(aes(xintercept = observed), colour = "red", linewidth = 0.6) +
  facet_wrap(~ feature, scales = "free", ncol = 5) +
  labs(x = "replicate value (red line: observed)", y = "count") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))
ggsave(sprintf("application-nhanes/figs/NHANES_ppcheck_%s.png", family), p,
       width = 11, height = 4.6, dpi = 300)
cat(sprintf("figs/NHANES_ppcheck_%s.png written\n", family))
