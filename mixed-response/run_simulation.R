## Fit the five copula models to the two simulated mixed-response data sets,
## four chains each, in parallel. Run from the repository root:
##
##   M=300000 CORES=8 Rscript mixed-response/run_simulation.R
##
## Optional environment variables:
##   M         iterations per chain (default 300000; the first third is burn-in)
##   CORES     parallel workers (default: all but one)
##   DATASETS  comma-separated subset of Gaussian,Gumbel
##   FAMILIES  comma-separated subset of Gaussian,Clayton,Frank,Gumbel,Independence
##   CHAINS    comma-separated chain indices (default 1,2,3,4)
##
## Each chain is written to output/mixed/fits/<data>_<family>_chain<k>.rds and
## a chain whose file already exists is skipped, so an interrupted batch can be
## resumed by rerunning the same command. Chain 1 starts at the initialisation
## of Section 3.1 of the paper; chains 2-4 start from that initialisation
## perturbed by N(0, 0.2^2) noise on each parameter and W drawn from its
## prior, so that the Gelman-Rubin diagnostic has dispersed starts to work with.

source("mixed-response/mixed_sampler.R")
suppressMessages(library(parallel))

M      <- as.integer(Sys.getenv("M", "300000"))
burn   <- floor(M / 3)
cores  <- as.integer(Sys.getenv("CORES", max(1L, detectCores() - 1L)))
split_ <- function(s, d) if (nzchar(s)) strsplit(s, ",")[[1]] else d
datasets <- split_(Sys.getenv("DATASETS"), c("Gaussian", "Gumbel"))
families <- split_(Sys.getenv("FAMILIES"), c("Gumbel", "Clayton", "Frank", "Gaussian", "Independence"))
chains   <- as.integer(split_(Sys.getenv("CHAINS"), 1:4))

dir.create("output/mixed/fits", recursive = TRUE, showWarnings = FALSE)
dir.create("output/mixed/logs", recursive = TRUE, showWarnings = FALSE)

jobs <- expand.grid(chain = chains, family = families, data = datasets,
                    stringsAsFactors = FALSE)
## Gumbel is the slowest family (its inverse h-function is numerical); it is
## listed first so that the long jobs start early
jobs <- jobs[order(match(jobs$family, c("Gumbel", "Clayton", "Frank", "Gaussian", "Independence"))), ]
jobs$file <- sprintf("output/mixed/fits/%s_%s_chain%d.rds", jobs$data, jobs$family, jobs$chain)
todo <- jobs[!file.exists(jobs$file), ]
cat(sprintf("%d of %d chains to run, M = %d, burn = %d, %d cores\n",
            nrow(todo), nrow(jobs), M, burn, cores))

one <- function(i) {
  j <- todo[i, ]
  log <- sprintf("output/mixed/logs/%s_%s_chain%d.log", j$data, j$family, j$chain)
  sink(log); on.exit(sink())
  d <- readRDS(sprintf("mixed-response/data/sim_%s.rds", j$data))
  t0 <- Sys.time()
  fit <- run_mixed(d$q, d$y, d$X, family = j$family, M = M, burn = burn,
                   seed = 1729 + 100 * j$chain, jitter = if (j$chain == 1) 0 else 0.2,
                   verbose = 25000)
  fit$data <- j$data; fit$chain <- j$chain
  fit$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(fit, j$file)
  cat(sprintf("done in %.1f min\n", fit$minutes))
  sprintf("%s %s chain %d: %.1f min", j$data, j$family, j$chain, fit$minutes)
}

if (nrow(todo) > 0) {
  res <- mclapply(seq_len(nrow(todo)), one, mc.cores = cores, mc.preschedule = FALSE)
  for (r in res) cat(if (inherits(r, "try-error")) paste("FAILED:", r) else r, "\n")
}
