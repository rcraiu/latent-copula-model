## Where does the sampler actually spend its time?
source("R/lvc_gamma_gauss.R")

NH <- readRDS("application-nhanes/nhanes_kidney.rds")
q <- NH$q; y <- NH$y; X <- NH$X

Rprof("/tmp/lvc.out", interval = 0.002, line.profiling = TRUE, memory.profiling = TRUE)
fit <- run_lvc(q, y, X, family = "Gaussian", M = 3000, burn = 1500,
               seed = 5, verbose = 0, thin_W = 0)
Rprof(NULL)

s <- summaryRprof("/tmp/lvc.out", memory = "none")
cat("=== total sampling time:", s$sampling.time, "sec for 3000 iterations\n\n")
cat("=== by self time (top 22) ===\n")
print(head(s$by.self, 22))
cat("\n=== by total time (top 14) ===\n")
print(head(s$by.total[s$by.total$total.pct < 99, ], 14))

## Count of primitive calls per iteration, for reference
n <- length(q)
cat("\nn =", n, "\n")
