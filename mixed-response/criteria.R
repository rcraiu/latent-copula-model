## Post-processing for the mixed-response simulation: posterior means,
## Gelman-Rubin diagnostics, and the model comparison criteria, both with W
## integrated out (the version the paper reports) and at the sampled W (the
## version originally defined). Run from the repository root after
## run_simulation.R:
##
##   CORES=8 Rscript mixed-response/criteria.R
##
## The W-marginalised criteria need one quadrature per retained draw and per
## observation, so the chains are thinned to DRAWS draws each (default 500,
## i.e. 2000 per model over four chains). Output goes to output/mixed/tables/
## as LaTeX fragments in the layout of Tables 2-7 of the paper, and to
## output/mixed/summary.rds.

source("mixed-response/mixed_model.R")
suppressMessages({ library(parallel); library(stableGR) })

cores <- as.integer(Sys.getenv("CORES", max(1L, detectCores() - 1L)))
ndraw <- as.integer(Sys.getenv("DRAWS", "500"))
datasets <- c("Gaussian", "Gumbel")
families <- c("Gaussian", "Clayton", "Frank", "Gumbel", "Independence")
dir.create("output/mixed/tables", recursive = TRUE, showWarnings = FALSE)

load_chains <- function(data, family) {
  files <- Sys.glob(sprintf("output/mixed/fits/%s_%s_chain*.rds", data, family))
  if (length(files) == 0) return(NULL)
  lapply(files, readRDS)
}

## ---- the criteria for one (data set, family) ------------------------------
one_model <- function(data, family) {
  fits <- load_chains(data, family)
  if (is.null(fits)) return(NULL)
  d <- readRDS(sprintf("mixed-response/data/sim_%s.rds", data))
  n <- length(d$q); pd <- ncol(d$X)
  gam <- fits[[1]]$gam
  st <- crit_state(n)
  for (f in fits) {
    idx <- round(seq(f$burn + 1, f$M, length.out = ndraw))
    for (t in idx) {
      th <- f$PAR[, t]
      beta <- th[seq_len(pd)]; p <- th[c("mud", "lamd", "muc", "lamc", "ssq")]
      xi <- c(th[["xi0"]], th[["xi1"]])
      L <- logdens_marg_W(d$q, d$y, as.numeric(d$X %*% beta), p, xi, family, gam)
      st <- crit_update(st, L)
    }
  }
  marg <- crit_waic(st); colnames(marg) <- c("joint", "d_g_c", "c_g_d")
  ## criteria at the sampled W, averaged over chains as in the original runs
  cond <- Reduce(`+`, lapply(fits, function(f) f$criteria_condW)) / length(fits)
  ## posterior means pooled over chains, and Gelman-Rubin over chains
  S <- lapply(fits, function(f) t(f$PAR[, (f$burn + 1):f$M]))
  keep <- if (family == "Independence") setdiff(colnames(S[[1]]), c("xi0", "xi1")) else colnames(S[[1]])
  S <- lapply(S, function(s) s[, keep, drop = FALSE])
  gr <- stable.GR(S, multivariate = TRUE); names(gr$psrf) <- colnames(S[[1]])
  list(data = data, family = family, marg = marg, cond = cond,
       means = colMeans(do.call(rbind, S)), psrf = gr$psrf, mpsrf = gr$mpsrf,
       draws = st$cnt, minutes = sapply(fits, function(f) f$minutes),
       acc = lapply(fits, function(f) f$acc))
}

grid <- expand.grid(family = families, data = datasets, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(grid)), function(i) one_model(grid$data[i], grid$family[i]),
                mc.cores = cores)
names(res) <- paste(grid$data, grid$family, sep = "_")
res <- Filter(Negate(is.null), res)
saveRDS(res, "output/mixed/summary.rds")

## ---- LaTeX tables ------------------------------------------------------------
fmt <- function(x, d = 4) ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
bold_min <- function(x, d = 2) {
  s <- fmt(x, d); s[which.min(x)] <- sprintf("\\textbf{%s}", s[which.min(x)]); s
}
parlab <- c(beta1 = "$\\beta_1$", beta2 = "$\\beta_2$", mud = "$\\mu_\\d$", lamd = "$\\lambda_\\d$",
            muc = "$\\mu_\\c$", lamc = "$\\lambda_\\c$", ssq = "$\\ssq$", xi0 = "$\\xi_0$", xi1 = "$\\xi_1$")
truth <- readRDS("mixed-response/data/sim_Gaussian.rds")$truth
truth_row <- c(truth$beta, truth$p, truth$xi)

for (data in datasets) {
  rr <- res[paste(data, families, sep = "_")]
  rr <- Filter(Negate(is.null), rr)
  if (length(rr) == 0) next
  fams <- sapply(rr, `[[`, "family")

  ## posterior means (Tables 2 and 5)
  M <- t(sapply(rr, function(r) r$means[names(parlab)][names(parlab)]))
  rownames(M) <- fams
  lines <- c("\\begin{tabular}{r|ccccccccc}", "  \\toprule",
             paste0(" & ", paste(parlab, collapse = " & "), " \\\\"), "  \\midrule",
             sprintf("  %s & %s \\\\", fams, apply(M, 1, function(v) paste(fmt(v), collapse = " & "))),
             "  \\midrule",
             sprintf("  True & %s \\\\", paste(fmt(truth_row, 1), collapse = " & ")),
             "   \\bottomrule", "\\end{tabular}")
  writeLines(lines, sprintf("output/mixed/tables/sim_%s_postmeans.tex", data))

  ## PSRFs (Tables 3 and 6)
  P <- t(sapply(rr, function(r) c(r$psrf[names(parlab)][names(parlab)], r$mpsrf)))
  lines <- c("\\begin{tabular}{r|cccccccccc}", "  \\toprule",
             paste0("& ", paste(c(parlab, "$\\bveta$"), collapse = " & "), " \\\\"), "  \\midrule",
             sprintf("  %s & %s \\\\", fams, apply(P, 1, function(v) paste(fmt(v), collapse = " & "))),
             "   \\bottomrule", "\\end{tabular}")
  writeLines(lines, sprintf("output/mixed/tables/sim_%s_psrf.tex", data))

  ## criteria (Tables 4 and 7): W integrated out, then at the sampled W
  crit_block <- function(get) {
    W  <- sapply(rr, function(r) get(r)["WAIC", "joint"])
    C1 <- sapply(rr, function(r) get(r)["WAIC", "d_g_c"])
    C2 <- sapply(rr, function(r) get(r)["WAIC", "c_g_d"])
    A  <- (C1 + C2) / 2
    c(sprintf("$\\text{WAIC}_\\text{est}$ & %s \\\\", paste(bold_min(W), collapse = " & ")),
      sprintf("$\\text{CWAIC}_{\\d \\mid \\c \\, \\text{est}}$ & %s \\\\", paste(bold_min(C1), collapse = " & ")),
      sprintf("$\\text{CWAIC}_{\\c \\mid \\d \\, \\text{est}}$ & %s \\\\", paste(bold_min(C2), collapse = " & ")),
      sprintf("Average $\\text{CWAIC}_{\\text{est}}$ & %s \\\\", paste(bold_min(A), collapse = " & ")))
  }
  lines <- c(sprintf("\\begin{tabular}{l|%s}", paste(rep("c", length(fams)), collapse = "")),
             "  \\toprule",
             paste0("& ", paste(fams, collapse = " & "), " \\\\"), "  \\midrule",
             sprintf("\\multicolumn{%d}{l}{\\emph{$W_i$ integrated out}} \\\\", length(fams) + 1),
             crit_block(function(r) r$marg), "  \\midrule",
             sprintf("\\multicolumn{%d}{l}{\\emph{evaluated at the sampled $W_i$}} \\\\", length(fams) + 1),
             crit_block(function(r) r$cond),
             "   \\bottomrule", "\\end{tabular}")
  writeLines(lines, sprintf("output/mixed/tables/sim_%s_criteria.tex", data))

  ## console summary
  cat(sprintf("\n=== data generated from the %s copula ===\n", data))
  show <- function(get, lab) {
    m <- sapply(rr, function(r) c(WAIC = get(r)["WAIC", "joint"], CWAIC_dc = get(r)["WAIC", "d_g_c"],
                                  CWAIC_cd = get(r)["WAIC", "c_g_d"],
                                  pwaic = get(r)["p_waic", "joint"]))
    colnames(m) <- fams
    cat(lab, "\n"); print(round(m, 2))
    cat("  selected by WAIC:", fams[which.min(m["WAIC", ])],
        " by CWAIC d|c:", fams[which.min(m["CWAIC_dc", ])],
        " by CWAIC c|d:", fams[which.min(m["CWAIC_cd", ])],
        " by average:", fams[which.min((m["CWAIC_dc", ] + m["CWAIC_cd", ]) / 2)], "\n")
  }
  show(function(r) r$marg, "-- W integrated out --")
  show(function(r) r$cond, "-- at the sampled W (as originally defined) --")
  cat("posterior means:\n"); print(round(M, 3))
  cat("max PSRF:", round(max(P[, 1:9], na.rm = TRUE), 4), " MPSRF:", round(P[, 10], 4), "\n")
}
