# =============================================================================
# 02_replicate.R
# Exact replication of Tables 1-3 in Ma & Cai (2025), "The Impact of Hospital
# Readmissions Reduction Program on Hospitalization and Home Care Use".
#
# Specification (as in the original code):
#   * sample: ages 60-70, people with >= 1 hospital discharge in the year
#   * local linear regression, uniform kernel, within the paper's bandwidths
#   * Y = a + b*1(age_c >= 0) + g*age_c + d*age_c*1(age_c >= 0) + year FE
#   * heteroskedasticity-robust (HC1) standard errors
# No methodological changes are made here; improvements are in 03 and 04.
# =============================================================================

d <- readRDS(ANALYSIS_FILE)
coef_map <- c("treated" = "Age 65+ (jump at cutoff)",
              "age_c" = "Age (months, centered)",
              "treated:age_c" = "Age 65+ x age")

# ---- Table 1: hospital discharges, 2013-2017 --------------------------------
t1 <- d %>% filter(post == 1, samp_orig_discharges)

bw_t1 <- rdrobust(y = t1$discharges, x = t1$age_c, c = 0)
message(sprintf("Table 1: rdrobust MSE-optimal bandwidth on current data = %.1f (paper used %.1f)",
                bw_t1$bws[1, 1], BW_ORIGINAL["discharges_post"]))

t1_local <- t1 %>% filter(abs(age_c) <= BW_ORIGINAL["discharges_post"])
m1a <- lm(discharges ~ treated * age_c + factor(YEAR), data = t1_local)

t1_local_c <- t1_local %>% filter(has_covs, !is.na(homecare_orig))
m1b <- lm(discharges ~ treated * age_c + factor(female) + factor(povcat) +
            factor(educ_group) + factor(race_group) + factor(region) +
            factor(YEAR) + homecare_orig, data = t1_local_c)

save_lm_table(
  list("Without controls" = m1a, "With controls" = m1b),
  list("Without controls" = vcovHC(m1a, type = "HC1"),
       "With controls"    = vcovHC(m1b, type = "HC1")),
  coef_map, "table1_replication_discharges",
  "Replication of Table 1: Hospital discharges (2013-2017)")

# ---- Table 2: home care days, 2013-2017 -------------------------------------
t2 <- d %>% filter(post == 1, samp_orig_homecare)

bw_t2 <- rdrobust(y = t2$hh_days[abs(t2$age_c) <= 60],
                  x = t2$age_c[abs(t2$age_c) <= 60], c = 0)
message(sprintf("Table 2: rdrobust bandwidth on current data = %.1f (paper used %.1f)",
                bw_t2$bws[1, 1], BW_ORIGINAL["homecare_post"]))

# In the original code both columns use the same covariate-complete sample.
t2_local <- t2 %>% filter(abs(age_c) <= BW_ORIGINAL["homecare_post"], has_covs)
m2a <- lm(hh_days ~ treated * age_c + factor(YEAR), data = t2_local)
m2b <- lm(hh_days ~ treated * age_c + factor(female) + factor(povcat) +
            factor(educ_group) + factor(race_group) + factor(region) +
            factor(YEAR), data = t2_local)

save_lm_table(
  list("Without controls" = m2a, "With controls" = m2b),
  list("Without controls" = vcovHC(m2a, type = "HC1"),
       "With controls"    = vcovHC(m2b, type = "HC1")),
  coef_map, "table2_replication_homecare",
  "Replication of Table 2: Home care days (2013-2017)")

# ---- Table 3: placebo with pre-HRRP data, 2007-2011 -------------------------
t3a <- d %>% filter(post == 0, samp_orig_discharges)
bw_t3a <- rdrobust(y = t3a$discharges, x = t3a$age_c, c = 0)
m3a <- lm(discharges ~ treated * age_c + factor(YEAR),
          data = t3a %>% filter(abs(age_c) <= BW_ORIGINAL["discharges_pre"]))

t3b <- d %>% filter(post == 0, samp_orig_homecare)
bw_t3b <- rdrobust(y = t3b$hh_days, x = t3b$age_c, c = 0)
m3b <- lm(hh_days ~ treated * age_c + factor(YEAR),
          data = t3b %>% filter(abs(age_c) <= BW_ORIGINAL["homecare_pre"]))

message(sprintf("Table 3: rdrobust bandwidths on current data = %.1f / %.1f (paper used %.1f / %.1f)",
                bw_t3a$bws[1, 1], bw_t3b$bws[1, 1],
                BW_ORIGINAL["discharges_pre"], BW_ORIGINAL["homecare_pre"]))

save_lm_table(
  list("Hospital discharges" = m3a, "Home care days" = m3b),
  list("Hospital discharges" = vcovHC(m3a, type = "HC1"),
       "Home care days"      = vcovHC(m3b, type = "HC1")),
  coef_map, "table3_replication_placebo",
  "Replication of Table 3: Placebo with pre-HRRP data (2007-2011)")

# ---- Compare with the numbers printed in the paper --------------------------
paper <- tibble::tribble(
  ~table,    ~model,                     ~paper_est, ~paper_se, ~paper_n,
  "Table 1", "Discharges, no controls",     -0.338,     0.141,     516,
  "Table 1", "Discharges, controls",        -0.384,     0.140,     508,
  "Table 2", "Home care days, no controls", -3.746,     6.322,     561,
  "Table 2", "Home care days, controls",    -3.144,     6.411,     561,
  "Table 3", "Placebo discharges",           0.276,     0.164,     437,
  "Table 3", "Placebo home care days",      -4.352,    11.031,     373)

models <- list(m1a, m1b, m2a, m2b, m3a, m3b)
replicated <- dplyr::bind_rows(lapply(models, function(m) {
  ct <- coeftest(m, vcov. = vcovHC(m, type = "HC1"))
  tibble::tibble(rep_est = ct["treated", 1], rep_se = ct["treated", 2], rep_n = nobs(m))
}))
comparison <- dplyr::bind_cols(paper, replicated) %>%
  mutate(matches = round(rep_est, 3) == paper_est &
                   round(rep_se, 3) == paper_se & rep_n == paper_n)
write.csv(comparison, file.path(DIR_TABLES, "replication_check.csv"), row.names = FALSE)
print(comparison, width = Inf)
if (!all(comparison$matches)) {
  message("Some estimates differ from the paper. Common causes: a different IPUMS ",
          "extract (sample selections / years), or a newer MEPS data release.")
}
