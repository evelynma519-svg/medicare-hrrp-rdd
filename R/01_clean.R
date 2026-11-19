# =============================================================================
# 01_clean.R
# Read the IPUMS MEPS extract and build one analysis file for both periods.
# Output: data/derived/analysis.rds
#
# Cleaning choices from the original paper are kept as *flags* rather than as
# hard filters, so that later scripts can reproduce the paper exactly
# (02_replicate.R) and also relax those choices (04_extensions.R).
# =============================================================================

# ---- Read raw extract -------------------------------------------------------
ddi_files <- list.files(DIR_RAW, pattern = "\\.xml$", full.names = TRUE)
if (length(ddi_files) == 0) {
  stop("No IPUMS DDI (.xml) file found in data/raw/. See README > Data.", call. = FALSE)
}
ddi_path <- ddi_files[which.max(file.mtime(ddi_files))]
message("Reading ", ddi_path)
raw <- ipumsr::read_ipums_micro(ipumsr::read_ipums_ddi(ddi_path), verbose = FALSE)

# Drop IPUMS value labels so variables behave as plain numbers.
raw[] <- lapply(raw, function(v) { attributes(v) <- NULL; v })

# Home care days: IPUMS now calls the person-level annual total HHTOTDAYSYR
# (older extracts: HHDAYS). Do NOT use HHDAYSMO, which is event-level.
if (!"HHDAYS" %in% names(raw) && "HHTOTDAYSYR" %in% names(raw)) {
  raw$HHDAYS <- raw$HHTOTDAYSYR
}

required <- c("YEAR", "AGE", "BIRTHYR", "BIRTHMO", "ENDRFY", "ENDRFM",
              "HPTOTDIS", "HHDAYS", "SEX", "EDUC", "POVCAT", "REGIONMEPS",
              "RACEA", "HIMCARE", "HEALTH", "ERTOTVIS", "RXPRMEDSNO", "OBEXPTOT")
absent <- setdiff(required, names(raw))
if (length(absent) > 0) {
  absent[absent == "HHDAYS"] <- "HHTOTDAYSYR"
  stop("Variables missing from the extract: ", paste(absent, collapse = ", "),
       ". See README > Data for the full variable list.", call. = FALSE)
}
for (v in c("MEPSID", "PERWEIGHT")) {
  if (!v %in% names(raw)) {
    warning(v, " not in extract: ",
            if (v == "MEPSID") "standard errors will not be clustered by person."
            else "survey-weighted robustness check will be skipped.", call. = FALSE)
  }
}

# ---- Build analysis file ----------------------------------------------------
analysis <- raw %>%
  filter(YEAR %in% c(PRE_YEARS, POST_YEARS),
         AGE >= AGE_MIN, AGE <= AGE_MAX) %>%
  mutate(
    # Running variable: age in months at the end of the survey reference
    # period, centered at 65 (same definition as the original code).
    age_in_months = (ENDRFY - BIRTHYR) * 12 + (ENDRFM - BIRTHMO),
    age_c         = age_in_months - CUTOFF_MONTHS,
    treated       = as.integer(age_c >= 0),
    # The original code did not screen birth/reference months for IPUMS
    # missing codes; 02_replicate.R keeps that, 03/04 drop invalid ages.
    valid_age     = BIRTHMO %in% 1:12 & ENDRFM %in% 1:12,
    post          = as.integer(YEAR %in% POST_YEARS),
    period        = ifelse(post == 1, "post (2013-2017)", "pre (2007-2011)"),

    # Someone just to the right of the cutoff turned 65 during the reference
    # year, so only part of their annual outcome is measured while eligible
    # for Medicare. Months of eligibility within the year (capped at 12):
    months_eligible = pmin(pmax(age_c, 0), 12),

    # Outcomes. Values outside these ranges are IPUMS missing / top codes.
    discharges  = ifelse(HPTOTDIS >= 0 & HPTOTDIS < 100, HPTOTDIS, NA_real_),
    any_hosp    = as.integer(discharges > 0),
    hh_days     = ifelse(HHDAYS >= 0 & HHDAYS < 365, HHDAYS, NA_real_),
    homecare_any = as.integer(hh_days > 0),
    # Original paper's version: built from raw HHDAYS, so IPUMS missing codes
    # (if any) count as "used home care". Kept only for exact replication.
    homecare_orig = as.integer(HHDAYS > 0),

    # First stage: has Medicare coverage (HIMCARE: 1 = no, 2 = yes)
    medicare    = case_when(HIMCARE == 2 ~ 1L, HIMCARE == 1 ~ 0L, TRUE ~ NA_integer_),

    # Covariates (same groupings as the original code)
    female      = case_when(SEX == 2 ~ 1L, SEX == 1 ~ 0L, TRUE ~ NA_integer_),
    educ_group  = case_when(
      EDUC >= 100 & EDUC <= 116 ~ "Less than high school",
      EDUC >= 200 & EDUC <= 202 ~ "High school graduate",
      EDUC >= 300 & EDUC <= 301 ~ "Some college",
      EDUC >= 302 & EDUC <= 303 ~ "Associate degree",
      EDUC == 400               ~ "Bachelor's degree",
      EDUC >= 500 & EDUC <= 599 ~ "Graduate degree",   # IPUMS now uses 510/520/530 for graduate detail
      TRUE ~ NA_character_),
    race_group  = case_when(
      RACEA == 100                ~ "White",
      RACEA == 200                ~ "Black",
      RACEA >= 300 & RACEA <= 350 ~ "American Indian/Alaska Native",
      RACEA >= 400 & RACEA <= 434 ~ "Asian/Pacific Islander",
      RACEA >= 500 & RACEA <= 617 ~ "Other/Mixed",
      TRUE ~ NA_character_),
    povcat      = ifelse(POVCAT %in% 1:5, POVCAT, NA),
    region      = ifelse(REGIONMEPS %in% 1:4, REGIONMEPS, NA),

    # Predetermined covariates for balance tests
    college     = as.integer(educ_group %in% c("Bachelor's degree", "Graduate degree")),
    white       = as.integer(race_group == "White"),
    black       = as.integer(race_group == "Black"),
    poor_nearpoor = as.integer(povcat %in% 1:2),

    # Other health/utilization outcomes (used as secondary outcomes)
    health_status = ifelse(HEALTH %in% 1:5, HEALTH, NA),      # 1 = excellent, 5 = poor
    er_visits     = ifelse(ERTOTVIS >= 0, ERTOTVIS, NA),
    rx_count      = ifelse(RXPRMEDSNO >= 0, RXPRMEDSNO, NA),
    office_exp    = ifelse(OBEXPTOT >= 0, OBEXPTOT, NA),

    has_covs    = !is.na(female) & !is.na(educ_group) & !is.na(race_group) &
                  !is.na(povcat) & !is.na(region),

    # Sample flags that reproduce the original paper's restrictions
    samp_orig_discharges = !is.na(discharges) & discharges > 0,
    samp_orig_homecare   = !is.na(hh_days) & !is.na(discharges) & discharges > 0
  )

keep_vars <- c("YEAR", intersect(c("MEPSID", "PERWEIGHT"), names(analysis)),
               "AGE", "age_in_months", "age_c", "treated", "valid_age", "post", "period",
               "months_eligible", "discharges", "any_hosp", "hh_days",
               "homecare_any", "homecare_orig", "medicare", "female", "educ_group", "race_group",
               "povcat", "region", "college", "white", "black", "poor_nearpoor",
               "health_status", "er_visits", "rx_count", "office_exp",
               "has_covs", "samp_orig_discharges", "samp_orig_homecare")
analysis <- analysis[, keep_vars]

saveRDS(analysis, ANALYSIS_FILE)

# ---- Sample construction log ------------------------------------------------
sample_log <- analysis %>%
  group_by(period) %>%
  summarise(
    persons_aged_60_70           = n(),
    with_valid_discharges        = sum(!is.na(discharges)),
    with_any_hospitalization     = sum(samp_orig_discharges),
    orig_sample_with_covariates  = sum(samp_orig_discharges & has_covs),
    .groups = "drop")
write.csv(sample_log, file.path(DIR_TABLES, "sample_construction.csv"), row.names = FALSE)
print(sample_log)
message("Saved ", ANALYSIS_FILE, " (", nrow(analysis), " rows)")
