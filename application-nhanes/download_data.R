## Download the NHANES 1999-2002 files needed for the kidney-function analysis.
##
## Run from the repository root:
##     Rscript application-nhanes/download_data.R
##
## The files are public and are fetched from the CDC. They total roughly 40 MB
## and are written to application-nhanes/data/, which is not tracked by git.
## Nothing is re-downloaded if it is already present, so the script is safe to
## re-run.

destdir <- "application-nhanes/data"
dir.create(destdir, recursive = TRUE, showWarnings = FALSE)

base99 <- "https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000"
base01 <- "https://wwwn.cdc.gov/Nchs/Nhanes/2001-2002"

files <- rbind(
  # cystatin C, measured on the surplus sera subsample
  data.frame(name = "SSCYST_A", url = file.path(base99, "SSCYST_A.XPT")),
  data.frame(name = "SSCYST_B", url = file.path(base01, "SSCYST_B.XPT")),
  # serum creatinine: LBXSCR in 1999-2000 (LAB18), LBDSCR in 2001-2002 (L40_B)
  data.frame(name = "LAB18",    url = file.path(base99, "LAB18.XPT")),
  data.frame(name = "L40_B",    url = file.path(base01, "L40_B.XPT")),
  # demographics, body measures, blood pressure
  data.frame(name = "DEMO",     url = file.path(base99, "DEMO.XPT")),
  data.frame(name = "DEMO_B",   url = file.path(base01, "DEMO_B.XPT")),
  data.frame(name = "BMX",      url = file.path(base99, "BMX.XPT")),
  data.frame(name = "BMX_B",    url = file.path(base01, "BMX_B.XPT")),
  data.frame(name = "BPX",      url = file.path(base99, "BPX.XPT")),
  data.frame(name = "BPX_B",    url = file.path(base01, "BPX_B.XPT"))
)

for (i in seq_len(nrow(files))) {
  dest <- file.path(destdir, paste0(files$name[i], ".xpt"))
  if (file.exists(dest) && file.size(dest) > 0) {
    cat(sprintf("%-10s already present (%.1f MB)\n",
                files$name[i], file.size(dest) / 1048576))
    next
  }
  cat(sprintf("%-10s downloading ... ", files$name[i]))
  ok <- tryCatch({
    utils::download.file(files$url[i], dest, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(dest) || file.size(dest) == 0) {
    if (file.exists(dest)) unlink(dest)
    stop("failed to download ", files$url[i],
         "\nCheck the URL in a browser; the CDC occasionally reorganises paths.")
  }
  cat(sprintf("%.1f MB\n", file.size(dest) / 1048576))
}

cat("\nall files present in ", destdir,
    "\nnext: Rscript application-nhanes/build_data.R\n", sep = "")
