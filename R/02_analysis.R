# ============================================================
# 02_analysis.R
# Core causal inference estimation functions.
#
# Implements the four estimators compared in Blackwell & Glynn (2018):
#   1.  Pooled OLS
#   2.  Unit fixed effects (within estimator)
#   3.  First differences
#   4.  IPW marginal structural model (MSM) — the paper's preferred method
#
# Also:
#   - Propensity score estimation (pooled logistic regression)
#   - Stabilised IPW weight construction
#   - Cluster-robust standard errors (clustered by country)
# ============================================================

source(here::here("R", "01_data_prep.R"))

dir.create(here::here("output"), showWarnings = FALSE)

# ============================================================
# 1. Propensity Score Estimation
# ============================================================

# ---- 1a. Denominator model: Pr(D_it | D_{i,t-1}, X_it) ----
#   Full model conditioning on time-varying covariates
ps_denom_formula <- democracy ~ democracy_lag + lgdppc_lag + education + trade

ps_denom_model <- glm(
  ps_denom_formula,
  data   = panel,
  family = binomial(link = "logit")
)

message("Denominator propensity-score model:")
print(summary(ps_denom_model))

# Predicted probabilities from denominator model
panel$ps_denom <- predict(ps_denom_model, type = "response")

# ---- 1b. Numerator model: Pr(D_it | D_{i,t-1}) ----
#   Stabilised weights: only condition on lagged treatment (no X)
ps_numer_formula <- democracy ~ democracy_lag

ps_numer_model <- glm(
  ps_numer_formula,
  data   = panel,
  family = binomial(link = "logit")
)

# Predicted probabilities from numerator model
panel$ps_numer <- predict(ps_numer_model, type = "response")

# ---- 1c. Compute per-period propensity for the observed treatment ----
panel <- panel %>%
  dplyr::mutate(
    # Pr(D_it = observed | ...) for each model
    p_obs_denom = dplyr::if_else(democracy == 1, ps_denom, 1 - ps_denom),
    p_obs_numer = dplyr::if_else(democracy == 1, ps_numer, 1 - ps_numer)
  )

# ---- 1d. Construct cumulative (product) stabilised weights ----
#   w_it = prod_{s=1}^{t} [numer_s / denom_s]
#   For practical reasons (numerical stability) we work in log-space
#   then exponentiate, and we ALSO compute simple (non-cumulative)
#   period-specific weights as used in many applied implementations.
panel <- panel %>%
  dplyr::group_by(country) %>%
  dplyr::arrange(year, .by_group = TRUE) %>%
  dplyr::mutate(
    # Log of per-period weight ratio
    log_w_period = log(p_obs_numer) - log(p_obs_denom),
    # Cumulative (product) weight — full MSM weight
    ipw_cumulative = exp(cumsum(log_w_period)),
    # Period-specific (non-cumulative) stabilised weight — simpler, commonly used
    ipw_period = p_obs_numer / p_obs_denom
  ) %>%
  dplyr::ungroup()

# ---- 1e. Trim extreme weights at 99th percentile ----
#   Extreme weights inflate variance; trimming (capping) at the 99th pct is standard practice
trim_99 <- quantile(panel$ipw_cumulative, 0.99, na.rm = TRUE)
panel <- panel %>%
  dplyr::mutate(
    ipw_cumulative_trim = pmin(ipw_cumulative, trim_99),
    # Normalise so mean weight = 1 (helps numerical stability)
    ipw_norm = ipw_cumulative_trim / mean(ipw_cumulative_trim, na.rm = TRUE)
  )

message("\nIPW weight summary (cumulative, trimmed):")
print(summary(panel$ipw_cumulative_trim))

# ============================================================
# 2. Define Covariate Formulas Used Across All Models
# ============================================================

outcome_formula_controls <- lgdppc ~ democracy + lgdppc_lag + education + trade
outcome_formula_simple   <- lgdppc ~ democracy

# ============================================================
# 3. Convert to plm Panel Data Frame
# ============================================================

panel_plm <- plm::pdata.frame(panel, index = c("country", "year"))

# ============================================================
# 4. Fit the Four Estimators
# ============================================================

# Helper: extract clustered SE coefficient table (for lm objects)
coef_table_cluster_lm <- function(model, cluster_var) {
  vcov_cl <- sandwich::vcovCL(model, cluster = cluster_var)
  lmtest::coeftest(model, vcov = vcov_cl)
}

# ---- 4a. Pooled OLS (with controls) ----
message("\n--- Pooled OLS ---")
ols_model <- lm(outcome_formula_controls, data = panel)
ols_coef  <- coef_table_cluster_lm(ols_model, cluster_var = panel$country)
print(ols_coef)

# ---- 4b. Fixed Effects (within, two-way) ----
#   Implemented as OLS with unit + year dummies for reliable cluster-SE handling
message("\n--- Fixed Effects (Two-Way) ---")
fe_model <- lm(
  lgdppc ~ democracy + lgdppc_lag + education + trade +
    factor(country) + factor(year),
  data = panel
)
# cluster-robust SE by country, then extract democracy row only
vcov_fe_cl <- sandwich::vcovCL(fe_model, cluster = panel$country)
fe_coef    <- lmtest::coeftest(fe_model, vcov = vcov_fe_cl)
message("Fixed Effects (two-way, democracy coefficient):")
print(fe_coef["democracy", , drop = FALSE])

# ---- 4c. First Differences ----
#   Use plm for automatic first-differencing
message("\n--- First Differences ---")
fd_model <- plm::plm(
  outcome_formula_simple,
  data  = panel_plm,
  model = "fd"
)
vcov_fd_cl <- plm::vcovHC(fd_model, method = "arellano", cluster = "group")
fd_coef    <- lmtest::coeftest(fd_model, vcov = vcov_fd_cl)
print(fd_coef)

# ---- 4d. IPW / Marginal Structural Model (preferred estimator) ----
#   Weighted two-way fixed effects using stabilised, trimmed IPW weights.
#   Implemented via weighted lm with factor() dummies (avoids plm weight
#   alignment issues and makes cluster SE computation straightforward).
message("\n--- IPW Marginal Structural Model ---")

# Drop any rows with missing or non-positive weights before fitting
msm_data <- panel %>%
  dplyr::filter(!is.na(ipw_norm), ipw_norm > 0)

msm_model <- lm(
  lgdppc ~ democracy + factor(country) + factor(year),
  data    = msm_data,
  weights = ipw_norm
)
vcov_msm_cl <- sandwich::vcovCL(msm_model, cluster = msm_data$country)
msm_coef    <- lmtest::coeftest(msm_model, vcov = vcov_msm_cl)
message("IPW/MSM (democracy coefficient):")
print(msm_coef["democracy", , drop = FALSE])

# ============================================================
# 5. Collect Results into a Tidy Table
# ============================================================

extract_row <- function(coef_obj, param = "democracy", label) {
  # coeftest returns a matrix; use positional indexing to be robust
  # across lm/glm/plm (column order: Estimate, SE, stat, p-value)
  idx <- which(rownames(coef_obj) == param)
  if (length(idx) == 0) return(NULL)
  data.frame(
    Estimator = label,
    Estimate  = coef_obj[idx, 1],
    Std.Error = coef_obj[idx, 2],
    t.value   = coef_obj[idx, 3],
    p.value   = coef_obj[idx, 4],
    stringsAsFactors = FALSE
  )
}

results_table <- dplyr::bind_rows(
  extract_row(ols_coef, label = "Pooled OLS"),
  extract_row(fe_coef,  label = "Two-Way FE"),
  extract_row(fd_coef,  label = "First Differences"),
  extract_row(msm_coef, label = "IPW / MSM (preferred)")
)

message("\n=== Summary of Democracy Effect Estimates ===")
print(results_table, digits = 4)

# ============================================================
# 6. Balance Diagnostics
# ============================================================

# Standardised mean differences for treated vs. control
# Before and after IPW weighting — used for Figure 2 type plots

smd <- function(x, trt, wt = NULL) {
  if (is.null(wt)) wt <- rep(1, length(x))
  m1 <- stats::weighted.mean(x[trt == 1], wt[trt == 1], na.rm = TRUE)
  m0 <- stats::weighted.mean(x[trt == 0], wt[trt == 0], na.rm = TRUE)
  v1 <- stats::var(x[trt == 1], na.rm = TRUE)
  v0 <- stats::var(x[trt == 0], na.rm = TRUE)
  (m1 - m0) / sqrt((v1 + v0) / 2)
}

covariates_for_balance <- c("lgdppc_lag", "education", "trade")

balance_df <- dplyr::bind_rows(lapply(covariates_for_balance, function(v) {
  data.frame(
    covariate = v,
    smd_unweighted = smd(panel[[v]], panel$democracy),
    smd_weighted   = smd(panel[[v]], panel$democracy, panel$ipw_norm),
    stringsAsFactors = FALSE
  )
}))

message("\nBalance diagnostics (standardised mean differences):")
print(balance_df)

# ============================================================
# 7. Save key objects for use in 03_tables.R and 04_figures.R
# ============================================================

message("\nSaving analysis objects to output/analysis_results.RData ...")
save(
  panel,
  ols_model, ols_coef,
  fe_model,  fe_coef,
  fd_model,  fd_coef,
  msm_model, msm_coef,
  results_table,
  balance_df,
  ps_denom_model,
  file = here::here("output", "analysis_results.RData")
)
message("Done.")
