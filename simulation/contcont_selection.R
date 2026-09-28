## Copula selection in the continuous-continuous (gamma-Gaussian) model, for
## the simulation section of the paper. Run from the repository root:
##
##   M=200000 CHAINS=2 CORES=8 Rscript simulation/contcont_selection.R
##
## Two data sets are generated with the design of the mixed-response
## simulation (covariates Unif(-3, 1) x Unif(-1, 1) standardised, beta = (1, 0.5),
## xi = (-1, 3), theta = kappa(0.1 (xi_0 + xi_1 W))), one with a Gaussian and
## one with a Gumbel copula, and margins
##   Q_i | W_i ~ Gamma(alpha = 5, mean exp(-1 + W_i)),   Y_i | W_i ~ N(W_i, 1).
## The sample size is N = 2000 (the size used for the experiments of Section
## 2.1): at n = 500 the likelihood is nearly flat along the ridge trading alpha
## against xi, the maximum likelihood estimate of xi sits far from the truth,
## and the copula family is barely identified.
##
## The five copula models are fitted with the C++ sampler, CHAINS chains of M
## iterations each, and scored two ways: by the criteria as originally defined,
## which evaluate p(q_i, y_i | eta, W_i) at the sampled latent variable and are
## accumulated online by the sampler, and by the versions with W_i integrated
## out by quadrature, which are what the paper reports. Both margins carry a free
## dispersion parameter here, so this is the case in which the two versions are
## expected to disagree most (Section 2.1).
##
## Environment variables: N (default 2000), M (default 200000, burn-in the
## first half), CHAINS (default 2), DRAWS per chain for the quadrature (default
## 250), CORES. Output: output/contcont/fits/, output/contcont/tables/,
## output/contcont/summary.rds; the data sets go to simulation/data/.

source("R/lvc_gamma_gauss.R")
source("R/lvc_fast.R")
suppressMessages({ library(parallel); library(stableGR) })

N      <- as.integer(Sys.getenv("N", "2000"))
M      <- as.integer(Sys.getenv("M", "200000"));  burn <- floor(M / 2)
nchain <- as.integer(Sys.getenv("CHAINS", "2"))
ndraw  <- as.integer(Sys.getenv("DRAWS", "250"))
cores  <- as.integer(Sys.getenv("CORES", max(1L, detectCores() - 1L)))
gscale <- 0.1
families <- c("Gaussian", "Clayton", "Frank", "Gumbel", "Independence")
for (d in c("output/contcont/fits", "output/contcont/tables", "simulation/data"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- data ----------------------------------------------------------------
## The generating values, and the same values in the identified
## parametrisation of Section 2.1 (||beta|| = 1, free sigma_W^2), which is
## what the sampler estimates: dividing W by ||beta|| multiplies the loadings
## and xi_1 by ||beta|| and sets sigma_W^2 = 1 / ||beta||^2.
gen <- list(beta = c(1, 0.5), muq = -1, lamq = 1, alpha = 5, muy = 0, lamy = 1,
            sigma = 1, xi = c(-1, 3))
nb  <- sqrt(sum(gen$beta^2))
truth <- c(muq = gen$muq, lamq = gen$lamq * nb, alpha = gen$alpha, muy = gen$muy,
           lamy = gen$lamy * nb, sigma2 = gen$sigma^2, xi0 = gen$xi[1], xi1 = gen$xi[2] * nb,
           beta1 = gen$beta[1] / nb, beta2 = gen$beta[2] / nb, s2W = 1 / nb^2)

make_data <- function(family, seed) {
  set.seed(seed)
  X <- scale(cbind(runif(N, -3, 1), runif(N, -1, 1)))
  X <- matrix(as.numeric(X), N, 2)
  W <- as.numeric(X %*% gen$beta) + rnorm(N)
  th <- theta_of_w(W, gen$xi, family, gscale)
  ## Rosenblatt construction through VineCopula's inverse h-function
  vf <- c(Gaussian = 1, Clayton = 3, Frank = 5, Gumbel = 4)[[family]]
  u2 <- runif(N)
  u1 <- VineCopula::BiCopHinv2(runif(N), u2, family = vf, par = th)
  q <- qgamma(u1, shape = gen$alpha, scale = exp(gen$muq + gen$lamq * W) / gen$alpha)
  y <- qnorm(u2, gen$muy + gen$lamy * W, gen$sigma)
  list(q = q, y = y, X = X, W = W, theta = th, family = family, gen = gen, truth = truth)
}
datasets <- c("Gaussian", "Gumbel")
for (f in datasets) {
  fn <- sprintf("simulation/data/contcont_%s.rds", f)
  if (!file.exists(fn)) saveRDS(make_data(f, seed = 2026 + match(f, datasets)), fn)
  d <- readRDS(fn)
  cat(sprintf("%s data: n = %d, Kendall tau(q, y) = %.3f, theta in [%.2f, %.2f], max q = %.1f\n",
              f, length(d$q), cor(d$q, d$y, method = "kendall"), min(d$theta), max(d$theta), max(d$q)))
}

## ---- fits ------------------------------------------------------------------
jobs <- expand.grid(chain = seq_len(nchain), family = families, data = datasets,
                    stringsAsFactors = FALSE)
jobs$file <- sprintf("output/contcont/fits/%s_%s_chain%d.rds", jobs$data, jobs$family, jobs$chain)
todo <- jobs[!file.exists(jobs$file), ]
cat(sprintf("%d of %d chains to run, M = %d\n", nrow(todo), nrow(jobs), M))
fit_one <- function(i) {
  j <- todo[i, ]
  d <- readRDS(sprintf("simulation/data/contcont_%s.rds", j$data))
  thin_W <- max(1L, floor((M - burn) / ndraw))
  fit <- run_lvc_fast(d$q, d$y, d$X, family = j$family, M = M, burn = burn,
                      gscale = gscale, seed = 1729 + 100 * j$chain, verbose = 0,
                      thin_W = thin_W, blockwise = TRUE, cross_block = TRUE)
  fit$data <- j$data; fit$chain <- j$chain
  saveRDS(fit, j$file)
  sprintf("%s %s chain %d done", j$data, j$family, j$chain)
}
if (nrow(todo) > 0)
  invisible(mclapply(seq_len(nrow(todo)), fit_one, mc.cores = cores, mc.preschedule = FALSE))

## ---- W integrated out --------------------------------------------------------
## p(q_i, y_i | eta) = int f(q_i, y_i | w) N(w; beta'x_i, s2W) dw and
## p(q_i | eta) likewise, by Simpson's rule. The integrand in w is proportional
## to the full conditional of W_i, so each observation gets its own grid,
## centred at the posterior mean of W_i and extending 8 posterior standard
## deviations either side (K nodes); the posterior of W_i pooled over draws
## of eta is wider than the conditional for any one draw, so the grid covers
## the peak wherever the draw puts it. The margin of y is in closed form.
## The step is checked against a grid twice as fine on every fiftieth draw, and
## the number of nodes is doubled and the model recomputed if the check fails
## (a badly misspecified copula can drive alpha and theta to values at which
## the integrand is extremely narrow).
logdens_marg <- function(q, y, m, sW, p, family, Wc, Ws, K = 401) {
  n <- length(q); aj <- numeric(n); aq <- numeric(n)
  h <- 16 * Ws / (K - 1)                                  # per-observation step
  wts <- c(1, rep(c(4, 2), length.out = K - 2), 1) / 3
  for (k in seq_len(K)) {
    W  <- Wc - 8 * Ws + (k - 1) * h
    wk <- wts[k] * h * dnorm(W, m, sW)
    cq <- comp_q(q, W, p); cy <- comp_y(y, W, p); th <- comp_th(W, p, family, gscale)
    aj <- aj + wk * exp(logdens_from(cq, cy, th, family, "joint"))
    aq <- aq + wk * exp(cq$lf)
  }
  lj <- log(pmax(aj, 1e-300)); lq <- log(pmax(aq, 1e-300))
  ly <- dnorm(y, p[["muy"]] + p[["lamy"]] * m,
              sqrt(exp(2 * p[["lsigma"]]) + p[["lamy"]]^2 * sW^2), log = TRUE)
  cbind(joint = lj, q_g_y = lj - ly, y_g_q = lj - lq)
}
crit_state <- function(n, k = 3) list(mp = matrix(0, n, k), ml = matrix(0, n, k),
                                       vl = matrix(0, n, k), cnt = 0)
crit_update <- function(s, L) {
  s$cnt <- s$cnt + 1; s$mp <- s$mp + (exp(L) - s$mp) / s$cnt
  dl <- L - s$ml; s$ml <- s$ml + dl / s$cnt; s$vl <- s$vl + dl * (L - s$ml); s
}
crit_waic <- function(s) {
  lppd <- colSums(log(pmax(s$mp, 1e-300))); pw <- colSums(s$vl) / s$cnt
  rbind(lppd = lppd, p_waic = pw, WAIC = -2 * lppd + 2 * pw)
}

one_model <- function(data, family) {
  fits <- lapply(Sys.glob(sprintf("output/contcont/fits/%s_%s_chain*.rds", data, family)), readRDS)
  d <- readRDS(sprintf("simulation/data/contcont_%s.rds", data)); n <- length(d$q)
  Wall <- do.call(cbind, lapply(fits, `[[`, "Wsamp"))
  Wc <- rowMeans(Wall); Ws <- pmax(apply(Wall, 1, sd), 0.05)
  K <- 401
  repeat {
    st <- crit_state(n); chk <- 0
    for (f in fits) {
      nd <- ncol(f$Wsamp); thin <- (f$M - f$burn) / nd
      idx <- f$burn + seq_len(nd) * thin
      for (j in seq_len(nd)) {
        th <- f$PAR[, idx[j]]; p <- th[1:8]
        m <- as.numeric(d$X %*% th[c("beta1", "beta2")]); sW <- sqrt(th[["s2W"]])
        L <- logdens_marg(d$q, d$y, m, sW, p, family, Wc, Ws, K = K)
        if (j %% 50 == 1) {                                 # quadrature check
          R <- logdens_marg(d$q, d$y, m, sW, p, family, Wc, Ws, K = 2 * K - 1)
          chk <- max(chk, max(abs(L - R)))
        }
        st <- crit_update(st, L)
      }
    }
    if (chk < 1e-3 || K >= 3201) break
    K <- 2 * K - 1
  }
  marg <- crit_waic(st)
  cond <- Reduce(`+`, lapply(fits, function(f) f$criteria[c("WAIC", "CWAIC_qy", "CWAIC_yq")])) / length(fits)
  S <- lapply(fits, function(f) t(f$PAR[, (f$burn + 1):f$M]))
  S <- lapply(S, function(s) {
    s <- cbind(s, alpha = exp(s[, "lalpha"]), sigma2 = exp(2 * s[, "lsigma"]))
    s[, setdiff(colnames(s), c("lalpha", "lsigma", if (family == "Independence") c("xi0", "xi1")))]
  })
  gr <- stable.GR(S, multivariate = TRUE); names(gr$psrf) <- colnames(S[[1]])
  list(data = data, family = family, marg = marg, cond = cond, quad_check = chk, K = K,
       means = colMeans(do.call(rbind, S)), psrf = gr$psrf, mpsrf = gr$mpsrf, draws = st$cnt)
}
grid <- expand.grid(family = families, data = datasets, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(grid)), function(i) one_model(grid$data[i], grid$family[i]),
                mc.cores = cores)
names(res) <- paste(grid$data, grid$family, sep = "_")
saveRDS(res, "output/contcont/summary.rds")

## ---- tables -----------------------------------------------------------------
fmt <- function(x, d = 4) ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
bold_min <- function(x, d = 2) { s <- fmt(x, d); s[which.min(x)] <- sprintf("\\textbf{%s}", s[which.min(x)]); s }
parlab <- c(beta1 = "$\\beta_1$", beta2 = "$\\beta_2$", s2W = "$\\sigma_W^2$", muq = "$\\mu_\\q$",
            lamq = "$\\lambda_\\q$", alpha = "$\\alpha$", muy = "$\\mu_\\y$", lamy = "$\\lambda_\\y$",
            sigma2 = "$\\ssq$", xi0 = "$\\xi_0$", xi1 = "$\\xi_1$")
for (data in datasets) {
  rr <- res[paste(data, families, sep = "_")]
  fams <- sapply(rr, `[[`, "family")
  Mn <- t(sapply(rr, function(r) r$means[names(parlab)][names(parlab)]))
  lines <- c(sprintf("\\begin{tabular}{r|%s}", paste(rep("c", length(parlab)), collapse = "")), "  \\toprule",
             paste0(" & ", paste(parlab, collapse = " & "), " \\\\"), "  \\midrule",
             sprintf("  %s & %s \\\\", fams, apply(Mn, 1, function(v) paste(fmt(v, 3), collapse = " & "))),
             "  \\midrule",
             sprintf("  True & %s \\\\", paste(fmt(truth[names(parlab)], 3), collapse = " & ")),
             "   \\bottomrule", "\\end{tabular}")
  writeLines(lines, sprintf("output/contcont/tables/contcont_%s_postmeans.tex", data))

  P <- t(sapply(rr, function(r) c(r$psrf[names(parlab)][names(parlab)], r$mpsrf)))
  lines <- c(sprintf("\\begin{tabular}{r|%s}", paste(rep("c", length(parlab) + 1), collapse = "")), "  \\toprule",
             paste0("& ", paste(c(parlab, "$\\bveta$"), collapse = " & "), " \\\\"), "  \\midrule",
             sprintf("  %s & %s \\\\", fams, apply(P, 1, function(v) paste(fmt(v), collapse = " & "))),
             "   \\bottomrule", "\\end{tabular}")
  writeLines(lines, sprintf("output/contcont/tables/contcont_%s_psrf.tex", data))

  block <- function(W, C1, C2) c(
    sprintf("$\\text{WAIC}_\\text{est}$ & %s \\\\", paste(bold_min(W), collapse = " & ")),
    sprintf("$\\text{CWAIC}_{\\q \\mid \\y \\, \\text{est}}$ & %s \\\\", paste(bold_min(C1), collapse = " & ")),
    sprintf("$\\text{CWAIC}_{\\y \\mid \\q \\, \\text{est}}$ & %s \\\\", paste(bold_min(C2), collapse = " & ")),
    sprintf("Average $\\text{CWAIC}_{\\text{est}}$ & %s \\\\", paste(bold_min((C1 + C2) / 2), collapse = " & ")))
  mW  <- sapply(rr, function(r) r$marg["WAIC", "joint"]); mC1 <- sapply(rr, function(r) r$marg["WAIC", "q_g_y"])
  mC2 <- sapply(rr, function(r) r$marg["WAIC", "y_g_q"])
  cW  <- sapply(rr, function(r) r$cond[["WAIC"]]); cC1 <- sapply(rr, function(r) r$cond[["CWAIC_qy"]])
  cC2 <- sapply(rr, function(r) r$cond[["CWAIC_yq"]])
  lines <- c(sprintf("\\begin{tabular}{l|%s}", paste(rep("c", length(fams)), collapse = "")), "  \\toprule",
             paste0("& ", paste(fams, collapse = " & "), " \\\\"), "  \\midrule",
             sprintf("\\multicolumn{%d}{l}{\\emph{$W_i$ integrated out}} \\\\", length(fams) + 1),
             block(mW, mC1, mC2), "  \\midrule",
             sprintf("\\multicolumn{%d}{l}{\\emph{evaluated at the sampled $W_i$}} \\\\", length(fams) + 1),
             block(cW, cC1, cC2), "   \\bottomrule", "\\end{tabular}")
  writeLines(lines, sprintf("output/contcont/tables/contcont_%s_criteria.tex", data))

  cat(sprintf("\n=== continuous data from the %s copula ===\n", data))
  out <- rbind(WAIC_marg = mW, CWAIC_qy_marg = mC1, CWAIC_yq_marg = mC2,
               WAIC_cond = cW, CWAIC_qy_cond = cC1, CWAIC_yq_cond = cC2,
               alpha = Mn[, "alpha"], sigma2 = Mn[, "sigma2"], s2W = Mn[, "s2W"],
               xi0 = Mn[, "xi0"], xi1 = Mn[, "xi1"],
               pwaic_marg = sapply(rr, function(r) r$marg["p_waic", "joint"]))
  colnames(out) <- fams; print(round(out, 2))
  cat("quadrature check (max abs log-density error):", signif(max(sapply(rr, `[[`, "quad_check")), 2), "\n")
  cat("max PSRF:", round(max(P[, -ncol(P)], na.rm = TRUE), 4), " MPSRF:", round(P[, ncol(P)], 4), "\n")
}
