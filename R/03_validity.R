# =============================================================================
# 03_validity.R
# Checks on the regression discontinuity design itself.
#
#   1. First stage: does Medicare coverage actually jump at 65?
#   2. Density of the running variable (manipulation / sample selection),
#      in the full sample AND in the paper's "hospitalized only" sample.
#   3. Balance of predetermined covariates at the cutoff.
#   4. Binned scatter plots of the main outcomes.
#
# Uses the post-HRRP period (2013-2017) and drops invalid ages.
# =============================================================================

d <- readRDS(ANALYSIS_FILE) %>% filter(valid_age)
post_all  <- d %>% filter(post == 1, !is.na(discharges))   # includes people with 0 discharges
post_hosp <- post_all %>% filter(samp_orig_discharges)      # paper's sample

cl <- function(df) if ("MEPSID" %in% names(df)) df$MEPSID else NULL

# ---- Binned scatter helper --------------------------------------------------
binned_plot <- function(df, y, ylab, title, window = 36, binwidth = 3) {
  df <- df[df$age_c >= -window & df$age_c < window & !is.na(df[[y]]), ]
  df$bin <- floor(df$age_c / binwidth) * binwidth + binwidth / 2
  bins <- df %>% group_by(bin) %>%
    summarise(mean_y = mean(.data[[y]]), n = n(), .groups = "drop")
  ggplot() +
    geom_point(data = bins, aes(bin, mean_y, size = n), alpha = 0.6) +
    geom_smooth(data = df[df$age_c < 0, ], aes(age_c, .data[[y]]),
                method = "lm", formula = y ~ x, color = "#2c7fb8", fill = "#a6bddb") +
    geom_smooth(data = df[df$age_c >= 0, ], aes(age_c, .data[[y]]),
                method = "lm", formula = y ~ x, color = "#d95f0e", fill = "#fdd0a2") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_size_area(max_size = 3, guide = "none") +
    labs(title = title, x = "Age in months relative to 65", y = ylab) +
    theme_minimal(base_size = 11)
}

# ---- 1. First stage ---------------------------------------------------------
fs <- rdrobust(y = post_all$medicare, x = post_all$age_c, c = 0,
               cluster = cl(post_all), vce = "hc1")
first_stage <- tidy_rd(fs, "Medicare coverage (first stage)")
save_fig(binned_plot(post_all, "medicare", "Share with Medicare",
                     "First stage: Medicare coverage around age 65"),
         "fig_first_stage")

# ---- 2. Density of the running variable -------------------------------------
dens_all  <- rddensity(X = post_all$age_c,  c = 0)
dens_hosp <- rddensity(X = post_hosp$age_c, c = 0)
density_tests <- tibble::tibble(
  sample = c("All respondents aged 60-70", "Hospitalized only (paper's sample)"),
  n      = c(nrow(post_all), nrow(post_hosp)),
  p_value_density_test = c(dens_all$test$p_jk, dens_hosp$test$p_jk))
write.csv(density_tests, file.path(DIR_TABLES, "validity_density_tests.csv"), row.names = FALSE)

hist_df <- dplyr::bind_rows(
  post_all  %>% transmute(age_c, sample = "All respondents"),
  post_hosp %>% transmute(age_c, sample = "Hospitalized only")) %>%
  filter(abs(age_c) <= 60)
p_hist <- ggplot(hist_df, aes(age_c)) +
  geom_histogram(binwidth = 3, boundary = 0, fill = "grey70", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#d95f0e") +
  facet_wrap(~ sample, scales = "free_y") +
  labs(title = "Density of the running variable", x = "Age in months relative to 65",
       y = "Observations") +
  theme_minimal(base_size = 11)
save_fig(p_hist, "fig_density", width = 9)

# ---- 3. Covariate balance ---------------------------------------------------
balance_vars <- c(female = "Female", college = "Bachelor's degree or higher",
                  white = "White", black = "Black",
                  poor_nearpoor = "Poor / near-poor (POVCAT 1-2)")
balance <- list()
for (smp in c("all", "hosp")) {
  df <- if (smp == "all") post_all else post_hosp
  for (v in names(balance_vars)) {
    ok <- !is.na(df[[v]])
    fit <- rdrobust(y = df[[v]][ok], x = df$age_c[ok], c = 0,
                    cluster = cl(df[ok, ]), vce = "hc1")
    balance[[paste(smp, v)]] <- tidy_rd(fit, balance_vars[[v]]) %>%
      mutate(sample = ifelse(smp == "all", "All respondents", "Hospitalized only"),
             .before = 1)
  }
}
balance <- dplyr::bind_rows(balance)
write.csv(balance, file.path(DIR_TABLES, "validity_covariate_balance.csv"), row.names = FALSE)

# ---- 4. Outcome plots -------------------------------------------------------
save_fig(binned_plot(post_hosp, "discharges", "Discharges per year",
                     "Hospital discharges, hospitalized only (paper's sample), 2013-2017"),
         "fig_rd_discharges_paper_sample")
save_fig(binned_plot(post_all, "any_hosp", "Share with any hospital stay",
                     "Any hospitalization, all respondents, 2013-2017"),
         "fig_rd_any_hosp")
save_fig(binned_plot(post_all, "discharges", "Discharges per year",
                     "Hospital discharges incl. zeros, all respondents, 2013-2017"),
         "fig_rd_discharges_all")

write.csv(first_stage, file.path(DIR_TABLES, "validity_first_stage.csv"), row.names = FALSE)
print(first_stage, width = Inf)
print(density_tests)
print(balance, width = Inf)
