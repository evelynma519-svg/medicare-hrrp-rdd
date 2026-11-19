# Medicare Eligibility, HRRP, and Hospital Use: Replication and Extension

An independent replication and extension of a regression discontinuity (RD) analysis of the Hospital Readmissions Reduction Program (HRRP), using the Medical Expenditure Panel Survey (MEPS).

**Original paper:** Evelyn Ma and Yuchen Cai (2025), *The Impact of Hospital Readmissions Reduction Program on Hospitalization and Home Care Use*, course paper for EC148 (Health Economics Research), Tufts University. The original analysis was co-authored. This repository, including all code, was written independently by Evelyn Ma to (1) reproduce the original results from the raw data and (2) address design issues identified after the paper was submitted.

---

## Research question

HRRP (in effect since October 2012) financially penalizes hospitals for excess 30-day readmissions among **Medicare** patients. Because Medicare eligibility begins at age 65, the original paper uses age as an RD running variable and asks whether annual hospital discharges drop discontinuously at 65 in the post-HRRP years (2013–2017).

The original paper found a reduction of about 0.34 discharges per year at age 65 among people with at least one hospitalization, and no change in home health care days.

## What this repository adds

| Step | Script | Why it matters |
|---|---|---|
| Exact replication of Tables 1–3 | `R/02_replicate.R` | Shows the published numbers can be reproduced from raw data; prints a side-by-side check against the paper (`output/tables/replication_check.csv`). |
| First stage | `R/03_validity.R` | Confirms that Medicare coverage actually jumps at 65 in this sample. |
| Density and covariate balance | `R/03_validity.R` | Tests the RD assumptions formally (`rddensity`, `rdrobust` on predetermined covariates), separately for all respondents and for the paper's hospitalized-only sample. |
| Modern RD inference | `R/04_extensions.R` | Triangular kernel, MSE-optimal bandwidth, robust bias-corrected confidence intervals (Calonico, Cattaneo & Titiunik 2014), standard errors clustered by person, because MEPS respondents appear in two consecutive years. |
| Keep people with zero discharges | `R/04_extensions.R` | The paper's sample only includes people hospitalized at least once. That conditions on the outcome: if Medicare changes *who* is hospitalized, the comparison at 65 is no longer like-for-like. The extension estimates the extensive margin (any hospitalization) and total discharges including zeros. |
| **Difference-in-discontinuities** | `R/04_extensions.R` | The jump at 65 bundles *every* effect of becoming eligible for Medicare (e.g. Card, Dobkin & Maestas 2008), not only HRRP. The paper's own placebo found a jump of +0.276 before HRRP. Subtracting the pre-HRRP jump (2007–2011) from the post-HRRP jump (2013–2017) isolates the change that coincides with HRRP (Grembi, Nannicini & Troiano 2016). |
| Donut RD | `R/04_extensions.R` | Age is measured at the end of the year, so people who turned 65 during the year are only partly exposed. They are dropped (Barreca et al. 2011). |
| Bandwidth sensitivity, survey weights | `R/04_extensions.R` | Shows whether conclusions depend on the bandwidth or on MEPS person weights. |
| Bad control removed | `R/04_extensions.R` | The original specification controls for home care use, which is itself a potential outcome of the policy. The extensions use only predetermined covariates. |

## Repository structure

```
medicare-hrrp-rdd/
├── README.md
├── run_all.R              # runs everything in order
├── R/
│   ├── 00_setup.R         # packages, paths, parameters, helpers
│   ├── 01_clean.R         # IPUMS extract -> data/derived/analysis.rds
│   ├── 02_replicate.R     # exact replication of the paper's Tables 1-3
│   ├── 03_validity.R      # first stage, density, balance, RD plots
│   └── 04_extensions.R    # rdrobust, extensive margin, diff-in-disc, donut, bandwidths
├── data/
│   ├── raw/               # put the IPUMS extract here (not tracked)
│   └── derived/           # created by 01_clean.R (not tracked)
└── output/
    ├── tables/            # CSV (and LaTeX if modelsummary is installed)
    └── figures/           # PNG
```

## Data

The data come from **IPUMS MEPS** (https://meps.ipums.org), which is free with registration. IPUMS terms of use do not allow redistribution, so the raw data are not included in this repository.

1. Create an extract with these **samples**: 2007, 2008, 2009, 2010, 2011, 2013, 2014, 2015, 2016, 2017.
2. Include these **variables**:

   | Purpose | Variables |
   |---|---|
   | Running variable | `AGE`, `BIRTHYR`, `BIRTHMO`, `ENDRFY`, `ENDRFM` |
   | Outcomes | `HPTOTDIS` (hospital discharges), `HHTOTDAYSYR` (annual days with home health visits; called `HHDAYS` in older extracts). Do not use `HHDAYSMO`, which is event-level. |
   | First stage | `HIMCARE` (Medicare coverage) |
   | Covariates | `SEX`, `EDUC`, `RACEA`, `POVCAT`, `REGIONMEPS` |
   | Secondary outcomes | `HEALTH`, `ERTOTVIS`, `RXPRMEDSNO`, `OBEXPTOT` |
   | Inference | `MEPSID` (person identifier), `PERWEIGHT` (person weight) |

   `YEAR` is added automatically. If `MEPSID` or `PERWEIGHT` are unavailable or named differently in your extract, the code still runs: it falls back to heteroskedasticity-robust standard errors and skips the weighted check, with a warning.
3. Download the data file (`.dat.gz`) **and** the DDI codebook (`.xml`), and put both in `data/raw/`.

## How to run

Requires R (≥ 4.1). Install the packages once:

```r
install.packages(c("ipumsr", "dplyr", "ggplot2", "rdrobust", "rddensity",
                   "sandwich", "lmtest", "tibble", "modelsummary"))
```

Then, from the repository root:

```bash
Rscript run_all.R
```

All tables are written to `output/tables/` and all figures to `output/figures/`. The full pipeline runs in a few minutes.

## Results

*To be filled in after running on the full data:*

- Replication check (paper vs. reproduced): `replication_check.csv`
- First stage, density, balance: `validity_*.csv`, `fig_first_stage.png`, `fig_density.png`
- RD estimates under the improved specification: `ext_rd_estimates_post.csv`
- Difference-in-discontinuities: `ext_diff_in_disc.csv`
- Bandwidth sensitivity: `fig_bandwidth_sensitivity.png`

## Notes and limitations

- **Running variable.** Age in months is measured at the end of each survey reference period (`ENDRFY`/`ENDRFM`), and outcomes are annual totals. That is why the donut specification is reported.
- **Outcome measure.** MEPS records discharges, not readmissions, and does not identify all six HRRP target conditions. Estimates refer to overall hospital use, not 30-day readmissions specifically.
- **No geographic detail.** Public-use MEPS has no state identifier, so state-level heterogeneity in HRRP penalties cannot be studied.
- **Interpretation of the diff-in-disc.** It attributes the change in the age-65 jump between periods to HRRP. Any other Medicare-specific change between 2007–2011 and 2013–2017 (e.g. other Affordable Care Act provisions affecting Medicare) would also be captured.

## References

- Barreca, A. I., Guldi, M., Lindo, J. M., & Waddell, G. R. (2011). Saving babies? Revisiting the effect of very low birth weight classification. *Quarterly Journal of Economics*, 126(4), 2117–2123.
- Calonico, S., Cattaneo, M. D., & Titiunik, R. (2014). Robust nonparametric confidence intervals for regression-discontinuity designs. *Econometrica*, 82(6), 2295–2326.
- Card, D., Dobkin, C., & Maestas, N. (2008). The impact of nearly universal insurance coverage on health care utilization: Evidence from Medicare. *American Economic Review*, 98(5), 2242–2258.
- Cattaneo, M. D., Jansson, M., & Ma, X. (2020). Simple local polynomial density estimators. *Journal of the American Statistical Association*, 115(531), 1449–1455.
- Grembi, V., Nannicini, T., & Troiano, U. (2016). Do fiscal rules matter? *American Economic Journal: Applied Economics*, 8(3), 1–30.
- Data: IPUMS MEPS, University of Minnesota. Please cite the data as requested on the IPUMS MEPS website.

## Contact

Evelyn Ma · evelynma519@gmail.com · [LinkedIn](https://www.linkedin.com/in/-evelynma-)
