# =============================================================================
# run_all.R  --  reproduce everything from the raw IPUMS extract.
# Usage (from the repository root):  Rscript run_all.R
# =============================================================================
source(file.path("R", "00_setup.R"))
source(file.path("R", "01_clean.R"))       # -> data/derived/analysis.rds
source(file.path("R", "02_replicate.R"))   # -> exact replication of Tables 1-3
source(file.path("R", "03_validity.R"))    # -> first stage, density, balance, plots
source(file.path("R", "04_extensions.R"))  # -> rdrobust, diff-in-disc, donut, bandwidths
message("Done. Tables in output/tables/, figures in output/figures/.")
