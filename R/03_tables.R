# ============================================================
# 03_tables.R
# Produce and save all tables from the replication.
#
# Corresponds to tables in Blackwell & Glynn (2018):
#   Table 1: Summary Statistics
#   Table 2: Propensity Score Model (denominator)
#   Table 3: Main Estimates — effect of democracy on log GDP pc
# ============================================================

source(here::here("R", "02_analysis.R"))

dir.create(here::here("output"), showWarnings = FALSE)

# ============================================================
# Table 1: Summary Statistics
# ============================================================

sum_stat <- panel %>%
  dplyr::select(
    `Log GDP per capita`   = lgdppc,
    `Democracy (D=1)`      = democracy,
    `Education`            = education,
    `Trade Openness`       = trade,
    `Lagged Log GDP pc`    = lgdppc_lag
  ) %>%
  tidyr::pivot_longer(everything(), names_to = "Variable", values_to = "value") %>%
  dplyr::group_by(Variable) %>%
  dplyr::summarise(
    N    = dplyr::n(),
    Mean = mean(value, na.rm = TRUE),
    SD   = sd(value,   na.rm = TRUE),
    Min  = min(value,  na.rm = TRUE),
    Max  = max(value,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 3)))

message("\n=== Table 1: Summary Statistics ===")
print(sum_stat)

# Save as CSV and LaTeX
write.csv(sum_stat, here::here("output", "table1_summary_stats.csv"), row.names = FALSE)

table1_kable <- kableExtra::kable(
  sum_stat,
  format  = "latex",
  booktabs = TRUE,
  caption = "Summary Statistics",
  label   = "tab:sumstats",
  digits  = 3
) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"))

writeLines(as.character(table1_kable),
           here::here("output", "table1_summary_stats.tex"))

message("Table 1 saved.")

# ============================================================
# Table 2: Propensity Score / Denominator Model Coefficients
# ============================================================

ps_coef_df <- broom::tidy(ps_denom_model, conf.int = TRUE) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  dplyr::select(
    Term      = term,
    Estimate  = estimate,
    `Std. Error` = std.error,
    `z value` = statistic,
    `Pr(>|z|)` = p.value,
    `95% CI Lower` = conf.low,
    `95% CI Upper` = conf.high
  )

message("\n=== Table 2: Propensity Score Model ===")
print(ps_coef_df)

write.csv(ps_coef_df, here::here("output", "table2_ps_model.csv"), row.names = FALSE)

table2_kable <- kableExtra::kable(
  ps_coef_df,
  format   = "latex",
  booktabs = TRUE,
  caption  = "Propensity Score Model (Logistic Regression, Denominator)",
  label    = "tab:psmodel",
  digits   = 4
) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"))

writeLines(as.character(table2_kable),
           here::here("output", "table2_ps_model.tex"))

message("Table 2 saved.")

# ============================================================
# Table 3: Main Results — Effect of Democracy on Log GDP pc
# ============================================================

# Format p-values
fmt_pval <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "< 0.001",
    p < 0.01  ~ paste0(round(p, 3), "**"),
    p < 0.05  ~ paste0(round(p, 3), "*"),
    TRUE      ~ as.character(round(p, 3))
  )
}

results_formatted <- results_table %>%
  dplyr::mutate(
    `95% CI` = paste0("[", round(Estimate - 1.96 * Std.Error, 3),
                      ", ", round(Estimate + 1.96 * Std.Error, 3), "]"),
    p.value  = fmt_pval(p.value),
    Estimate = round(Estimate, 4),
    Std.Error = round(Std.Error, 4)
  ) %>%
  dplyr::select(Estimator, Estimate, `Std. Error` = Std.Error, `95% CI`, `p-value` = p.value)

message("\n=== Table 3: Effect of Democracy on Log GDP per Capita ===")
print(results_formatted)

write.csv(results_formatted, here::here("output", "table3_main_results.csv"), row.names = FALSE)

table3_kable <- kableExtra::kable(
  results_formatted,
  format   = "latex",
  booktabs = TRUE,
  caption  = paste0(
    "Effect of Democracy on Log GDP per Capita under Alternative Estimators.",
    " Cluster-robust standard errors (by country) in parentheses.",
    " * p < 0.05, ** p < 0.01, *** p < 0.001."
  ),
  label    = "tab:mainresults"
) %>%
  kableExtra::kable_styling(latex_options = c("hold_position")) %>%
  kableExtra::row_spec(nrow(results_formatted), bold = TRUE)  # highlight MSM row

writeLines(as.character(table3_kable),
           here::here("output", "table3_main_results.tex"))

message("Table 3 saved.")

# ============================================================
# Table 4: Balance Diagnostics
# ============================================================

balance_formatted <- balance_df %>%
  dplyr::rename(
    Covariate                    = covariate,
    `SMD (Unweighted)`           = smd_unweighted,
    `SMD (IPW Weighted)`         = smd_weighted
  ) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4)))

message("\n=== Table 4: Balance Diagnostics ===")
print(balance_formatted)

write.csv(balance_formatted, here::here("output", "table4_balance.csv"), row.names = FALSE)

table4_kable <- kableExtra::kable(
  balance_formatted,
  format   = "latex",
  booktabs = TRUE,
  caption  = paste0(
    "Standardised Mean Differences (SMD) Before and After IPW Weighting.",
    " Values close to 0 indicate good balance between treated and control units."
  ),
  label    = "tab:balance"
) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"))

writeLines(as.character(table4_kable),
           here::here("output", "table4_balance.tex"))

message("Table 4 saved.")
message("\nAll tables written to output/")
