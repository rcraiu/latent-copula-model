## Download the NHANES files needed for the kidney-function analysis.
##
## Run from the repository root:
##     Rscript application-nhanes/download_data.R
##
## The files are public and are fetched from the CDC. They total roughly 45 MB
## and are written to application-nhanes/data/, which is not tracked by git.
## Nothing is re-downloaded if it is already present, so the script is safe to
## re-run.
##
## Note on the URLs: the CDC serves these from
##     .../Nchs/Data/Nhanes/Public/<cycle begin year>/DataFiles/<name>.xpt
## and *not* from the .../Nchs/Nhanes/<year range>/<name>.XPT form that older
## documentation and many papers cite; that form now returns an HTML "page not
## found" document with status 200. Each download is therefore checked for the
## SAS transport file signature rather than trusting the status code, since a
## silently saved HTML error page is the failure mode to guard against.

destdir <- "application-nhanes/data"
dir.create(destdir, recursive = TRUE, showWarnings = FALSE)

base <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"

## The surplus sera file spans 1999-2004, but the analysis uses the 1999-2000
## and 2001-2002 cycles, so only those covariate files are needed.
files <- rbind(
  # Cystatin C and beta-2 microglobulin, 1999-2004 surplus sera. This single
  # file supplies both surrogate markers.
  data.frame(name = "SSCARD_A", year = 1999),
  # demographics, body measures, blood pressure, by cycle
  data.frame(name = "DEMO",   year = 1999),
  data.frame(name = "DEMO_B", year = 2001),
  data.frame(name = "BMX",    year = 1999),
  data.frame(name = "BMX_B",  year = 2001),
  data.frame(name = "BPX",    year = 1999),
  data.frame(name = "BPX_B",  year = 2001),
  # serum creatinine, used only to stage kidney disease for the exclusion
  data.frame(name = "LAB18",  year = 1999),
  data.frame(name = "L40_B",  year = 2001)
)

## A SAS transport file begins with this header record.
is_xport <- function(path) {
  if (!file.exists(path) || file.size(path) < 80) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con))
  grepl("HEADER RECORD", rawToChar(readBin(con, "raw", 80)), fixed = TRUE)
}

for (i in seq_len(nrow(files))) {
  nm   <- files$name[i]
  dest <- file.path(destdir, paste0(nm, ".xpt"))
  if (is_xport(dest)) {
    cat(sprintf("%-10s already present (%.1f MB)\n", nm, file.size(dest) / 1048576))
    next
  }
  url <- sprintf("%s/%d/DataFiles/%s.xpt", base, files$year[i], nm)
  cat(sprintf("%-10s downloading ... ", nm))
  try(utils::download.file(url, dest, mode = "wb", quiet = TRUE), silent = TRUE)
  if (!is_xport(dest)) {
    if (file.exists(dest)) unlink(dest)          # do not leave an error page
    stop("failed to download a valid transport file from\n  ", url,
         "\nOpen that URL in a browser to check whether the CDC has moved it.")
  }
  cat(sprintf("%.1f MB\n", file.size(dest) / 1048576))
}

cat("\nall files present and valid in ", destdir,
    "\nnext: Rscript application-nhanes/build_data.R\n", sep = "")
