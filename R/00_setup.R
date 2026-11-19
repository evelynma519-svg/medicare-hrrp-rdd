# =============================================================================
# 00_setup.R
# Packages, paths, analysis parameters, and small helper functions.
# Sourced by every other script (via run_all.R). Run everything from the
# repository root so the relative paths below resolve.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
# Packages are NOT installed automatically. If any are missing, the script
# stops and prints the exact install.packages() call to run once.
PKGS <- c("ipumsr", "dplyr", "ggplot2", "rdrobust", "rddensity",
          "sandwich", "lmtest", "tibble")
PKGS_OPTIONAL <- c("modelsummary")   # only used for LaTeX tables

missing_pkgs <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Run once:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))", call. = FALSE)
}
for (p in PKGS) suppressPackageStartupMessages(library(p, character.only = TRUE))
HAS_MODELSUMMARY <- requireNamespace("modelsummary", quietly = TRUE)

# ---- Paths ------------------------------------------------------------------
if (!file.exists("README.md") || !dir.exists("R")) {
  stop("Please run the scripts from the repository root (the folder with README.md).",
       call. = FALSE)
}
DIR_RAW     <- file.path("data", "raw")
DIR_DERIVED <- file.path("data", "derived")
DIR_TABLES  <- file.path("output", "tables")
DIR_FIGS    <- file.path("output", "figures")
for (d in c(DIR_DERIVED, DIR_TABLES, DIR_FIGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
ANALYSIS_FILE <- file.path(DIR_DERIVED, "analysis.rds")

# ---- Analysis parameters ----------------------------------------------------
CUTOFF_MONTHS <- 65 * 12          # age 65 = 780 months
AGE_MIN       <- 60               # sample window in years (as in the paper)
AGE_MAX       <- 70
POST_YEARS    <- 2013:2017        # HRRP in effect (began Oct 2012)
PRE_YEARS     <- 2007:2011        # pre-HRRP placebo period

# Bandwidths reported in the original paper (Ma & Cai 2025). They were chosen
# by rdrobust and then typed in by hand; 02_replicate.R uses them so the
# original tables can be reproduced exactly, and also prints what rdrobust
# selects on the current data so any discrepancy is visible.
BW_ORIGINAL <- c(discharges_post = 18.8, homecare_post = 20.6,
                 discharges_pre  = 17.5, homecare_pre  = 15.1)

set.seed(20250506)

# ---- Helpers ----------------------------------------------------------------

# Tidy one rdrobust fit into a one-row data frame. Reports the conventional
# point estimate with its SE, plus the robust bias-corrected CI and p-value
# (Calonico, Cattaneo & Titiunik 2014), which is what should be used for
# inference with an MSE-optimal bandwidth.
tidy_rd <- function(fit, label) {
  tibble::tibble(
    spec       = label,
    estimate   = unname(fit$coef[1, 1]),
    se_conv    = unname(fit$se[1, 1]),
    ci_rb_low  = unname(fit$ci[3, 1]),
    ci_rb_high = unname(fit$ci[3, 2]),
    p_rb       = unname(fit$pv[3, 1]),
    bw_h       = unname(fit$bws[1, 1]),
    n_left     = unname(fit$N_h[1]),
    n_right    = unname(fit$N_h[2])
  )
}

# Tidy one coefficient from an lm fit with a given variance matrix.
tidy_lm_coef <- function(fit, vcov, term, label) {
  ct <- lmtest::coeftest(fit, vcov. = vcov)
  tibble::tibble(
    spec     = label,
    term     = term,
    estimate = unname(ct[term, 1]),
    se       = unname(ct[term, 2]),
    p        = unname(ct[term, 4]),
    n        = stats::nobs(fit)
  )
}

# Triangular kernel weights (the rdrobust default) for a bandwidth h.
tri_weights <- function(x, h) pmax(0, 1 - abs(x) / h)

# Cluster-robust variance if a person identifier is available, else HC1.
vcov_cluster <- function(fit, data) {
  if ("MEPSID" %in% names(data) && !all(is.na(data$MEPSID))) {
    sandwich::vcovCL(fit, cluster = data$MEPSID, type = "HC1")
  } else {
    sandwich::vcovHC(fit, type = "HC1")
  }
}

# Save a list of lm models as a regression table: always a CSV of the
# requested coefficients; additionally a LaTeX table if modelsummary is
# installed.
save_lm_table <- function(models, vcovs, coef_map, file_stub, title,
                          add_rows = NULL) {
  rows <- list()
  for (nm in names(models)) {
    ct <- lmtest::coeftest(models[[nm]], vcov. = vcovs[[nm]])
    keep <- intersect(names(coef_map), rownames(ct))
    rows[[nm]] <- tibble::tibble(model = nm, term = keep,
                                 estimate = ct[keep, 1], se = ct[keep, 2],
                                 p = ct[keep, 4], n = stats::nobs(models[[nm]]))
  }
  out <- dplyr::bind_rows(rows)
  utils::write.csv(out, file.path(DIR_TABLES, paste0(file_stub, ".csv")), row.names = FALSE)
  if (HAS_MODELSUMMARY) {
    modelsummary::modelsummary(
      models, coef_map = coef_map, vcov = unname(vcovs),
      stars = c("*" = .05, "**" = .01, "***" = .001),
      gof_omit = "IC|Log|Adj|F|RMSE|Std.Errors",
      add_rows = add_rows, title = title, escape = FALSE,
      output = file.path(DIR_TABLES, paste0(file_stub, ".tex"))
    )
  }
  invisible(out)
}

save_fig <- function(plot, file_stub, width = 7, height = 4.5) {
  ggplot2::ggsave(file.path(DIR_FIGS, paste0(file_stub, ".png")), plot,
                  width = width, height = height, dpi = 300)
}
