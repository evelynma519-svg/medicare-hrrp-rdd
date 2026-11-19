# =============================================================================
# 04_extensions.R
# Extensions to the original analysis.
#
#   E1. Modern RD inference: rdrobust estimates with a triangular kernel,
#       MSE-optimal bandwidth, robust bias-corrected CIs, SEs clustered by
#       person (MEPS respondents appear in two consecutive years), and
#       covariates entered without the post-treatment "home care" control.
#   E2. Extensive margin: keep people with zero discharges (the paper drops
#       them, which conditions the sample on the outcome).
#   E3. Difference-in-discontinuities: the age-65 jump bundles every effect
#       of becoming Medicare-eligible. Subtracting the pre-HRRP jump
#       (2007-2011) from the post-HRRP jump (2013-2017) isolates the change
#       in the jump that coincides with HRRP.
#   E4. Donut RD: people just right of the cutoff turned 65 during the year,
#       so their annual outcome is only partly "treated". Drop them.
#   E5. Bandwidth sensitivity for the RD and diff-in-disc estimates.
#   E6. MEPS person weights.
# =============================================================================

d <- readRDS(ANALYSIS_FILE) %>% filter(valid_age, !is.na(discharges))
has_id <- "MEPSID" %in% names(d)
has_wt <- "PERWEIGHT" %in% names(d)
cl <- function(df) if (has_id) df$MEPSID else NULL

covs_matrix <- function(df) {
  model.matrix(~ factor(female) + factor(povcat) + factor(educ_group) +
                 factor(race_group) + factor(region) + factor(YEAR), data = df)[, -1, drop = FALSE]
}

rd <- function(df, y, label, covs = FALSE, weights = NULL, h = NULL) {
  df <- df[!is.na(df[[y]]), ]
  if (covs) df <- df[df$has_covs, ]
  fit <- rdrobust(y = df[[y]], x = df$age_c, c = 0,
                  covs = if (covs) covs_matrix(df) else NULL,
                  cluster = cl(df), vce = "hc1",
                  weights = if (!is.null(weights)) df[[weights]] else NULL,
                  h = h)
  tidy_rd(fit, label)
}

post      <- d %>% filter(post == 1)
DONUT_MONTHS <- 12
donut <- function(df) df %>% filter(!(age_c >= 0 & age_c < DONUT_MONTHS))
post_hosp <- post %>% filter(samp_orig_discharges)

# ---- E1, E2, E4, E6: RD estimates for the post-HRRP period ------------------
rd_rows <- list(
  rd(post_hosp, "discharges", "Discharges | hospitalized (paper's sample)"),
  rd(post_hosp, "discharges", "Discharges | hospitalized, with covariates", covs = TRUE),
  rd(post,      "discharges", "Discharges, incl. zeros (all respondents)"),
  rd(post,      "discharges", "Discharges, incl. zeros, with covariates", covs = TRUE),
  rd(post,      "any_hosp",   "Any hospitalization (extensive margin)"),
  rd(post,      "any_hosp",   "Any hospitalization, with covariates", covs = TRUE),
  # Donut RD: the hole (0-11 months after 65) is wider than a data-driven
  # bandwidth can be, so the bandwidth is fixed at 24 and 36 months.
  rd(donut(post), "discharges",
     "Donut (h = 24): discharges incl. zeros", h = 24),
  rd(donut(post), "discharges",
     "Donut (h = 36): discharges incl. zeros", h = 36),
  rd(donut(post_hosp), "discharges",
     "Donut (h = 24): discharges | hospitalized", h = 24),
  rd(donut(post_hosp), "discharges",
     "Donut (h = 36): discharges | hospitalized", h = 36),
  rd(post,      "hh_days",    "Home care days, all respondents"))
if (has_wt) {
  rd_rows <- c(rd_rows, list(
    rd(post_hosp, "discharges", "Discharges | hospitalized, weighted", weights = "PERWEIGHT"),
    rd(post,      "discharges", "Discharges, incl. zeros, weighted",   weights = "PERWEIGHT")))
}
rd_table <- dplyr::bind_rows(rd_rows)
write.csv(rd_table, file.path(DIR_TABLES, "ext_rd_estimates_post.csv"), row.names = FALSE)
print(rd_table, width = Inf)

# ---- E3: Difference-in-discontinuities --------------------------------------
# Pooled local linear regression with triangular weights:
#   Y = a + b1*T + g1*age + d1*T*age
#         + b2*T*Post + g2*age*Post + d2*T*age*Post + year FE + e
# b1 = jump at 65 before HRRP, b1 + b2 = jump after HRRP, b2 = diff-in-disc.
# (Post itself is absorbed by the year fixed effects.)
diff_in_disc <- function(df, y, h, label) {
  df <- df[!is.na(df[[y]]) & abs(df$age_c) <= h, ]
  df$w <- tri_weights(df$age_c, h)
  f <- as.formula(paste(y, "~ treated * age_c + treated:post + age_c:post +",
                        "treated:age_c:post + factor(YEAR)"))
  fit <- lm(f, data = df, weights = w)
  V <- vcov_cluster(fit, df)
  b <- coef(fit)
  post_jump    <- b["treated"] + b["treated:post"]
  post_jump_se <- sqrt(V["treated", "treated"] + V["treated:post", "treated:post"] +
                       2 * V["treated", "treated:post"])
  dd_se <- sqrt(V["treated:post", "treated:post"])
  tibble::tibble(
    spec = label, bandwidth = h, n = nobs(fit),
    jump_pre  = unname(b["treated"]),   se_pre  = sqrt(V["treated", "treated"]),
    jump_post = unname(post_jump),      se_post = unname(post_jump_se),
    diff_in_disc = unname(b["treated:post"]), se_did = dd_se,
    p_did = 2 * pnorm(-abs(unname(b["treated:post"]) / dd_se)))
}

# Bandwidth: MSE-optimal bandwidth for the post-period RD of the same outcome.
mse_bw <- function(df, y) {
  df <- df[!is.na(df[[y]]), ]
  rdrobust(y = df[[y]], x = df$age_c, c = 0, cluster = cl(df), vce = "hc1")$bws[1, 1]
}
hosp <- d %>% filter(samp_orig_discharges)
did_table <- dplyr::bind_rows(
  diff_in_disc(hosp, "discharges", mse_bw(post_hosp, "discharges"),
               "Discharges | hospitalized (paper's sample)"),
  diff_in_disc(d, "discharges", mse_bw(post, "discharges"),
               "Discharges, incl. zeros"),
  diff_in_disc(d, "any_hosp", mse_bw(post, "any_hosp"),
               "Any hospitalization"))
write.csv(did_table, file.path(DIR_TABLES, "ext_diff_in_disc.csv"), row.names = FALSE)
print(did_table, width = Inf)

# ---- E5: Bandwidth sensitivity ----------------------------------------------
bw_grid <- seq(6, 48, by = 3)
local_rd <- function(df, y, h) {
  df <- df[!is.na(df[[y]]) & abs(df$age_c) <= h, ]
  df$w <- tri_weights(df$age_c, h)
  fit <- lm(as.formula(paste(y, "~ treated * age_c + factor(YEAR)")), data = df, weights = w)
  V <- vcov_cluster(fit, df)
  c(est = unname(coef(fit)["treated"]), se = sqrt(V["treated", "treated"]))
}
sens <- dplyr::bind_rows(lapply(bw_grid, function(h) {
  rd_post <- local_rd(post_hosp, "discharges", h)
  dd      <- diff_in_disc(hosp, "discharges", h, "")
  tibble::tibble(
    bandwidth = h,
    estimator = c("RD, post-HRRP (paper's sample)", "Diff-in-disc (paper's sample)"),
    estimate  = c(rd_post[["est"]], dd$diff_in_disc),
    se        = c(rd_post[["se"]],  dd$se_did))
}))
write.csv(sens, file.path(DIR_TABLES, "ext_bandwidth_sensitivity.csv"), row.names = FALSE)

p_sens <- ggplot(sens, aes(bandwidth, estimate)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_ribbon(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se),
              fill = "grey85") +
  geom_line() + geom_point(size = 1.5) +
  facet_wrap(~ estimator) +
  labs(title = "Sensitivity of the estimated jump at 65 to the bandwidth",
       subtitle = "Local linear, triangular kernel; shaded band = 95% CI",
       x = "Bandwidth (months on each side of 65)", y = "Estimated effect on discharges") +
  theme_minimal(base_size = 11)
save_fig(p_sens, "fig_bandwidth_sensitivity", width = 9)
