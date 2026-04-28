# ============================================================
# 04_figures.R
# Produce and save all figures from the replication.
#
# Corresponds to figures in Blackwell & Glynn (2018):
#   Figure 1: Distribution of IPW weights
#   Figure 2: Covariate balance before and after weighting (love plot)
#   Figure 3: Democracy effect estimates across estimators
#   Figure 4: Propensity scores for treated and control units
# ============================================================

source(here::here("R", "02_analysis.R"))

dir.create(here::here("output"), showWarnings = FALSE)

theme_replication <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90", colour = NA),
      legend.position  = "bottom"
    )
}

# ============================================================
# Figure 1: Distribution of IPW Weights
# ============================================================

fig1 <- ggplot2::ggplot(panel, ggplot2::aes(x = ipw_cumulative_trim)) +
  ggplot2::geom_histogram(
    bins  = 50,
    fill  = "steelblue",
    color = "white",
    alpha = 0.8
  ) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "red", linewidth = 0.7) +
  ggplot2::labs(
    x       = "Stabilised IPW Weight (trimmed at 99th percentile)",
    y       = "Count",
    title   = "Figure 1: Distribution of Inverse Probability Weights",
    caption = "Red dashed line at weight = 1 (no reweighting)."
  ) +
  ggplot2::scale_x_continuous(labels = scales::comma) +
  theme_replication()

ggplot2::ggsave(
  here::here("output", "figure1_ipw_weights.pdf"),
  fig1, width = 7, height = 4.5
)
ggplot2::ggsave(
  here::here("output", "figure1_ipw_weights.png"),
  fig1, width = 7, height = 4.5, dpi = 150
)
message("Figure 1 saved.")

# ============================================================
# Figure 2: Love Plot — Covariate Balance
# ============================================================

balance_long <- balance_df %>%
  tidyr::pivot_longer(
    cols      = c(smd_unweighted, smd_weighted),
    names_to  = "Type",
    values_to = "SMD"
  ) %>%
  dplyr::mutate(
    Type = dplyr::recode(Type,
      smd_unweighted = "Unweighted",
      smd_weighted   = "IPW Weighted"
    ),
    Type      = factor(Type, levels = c("Unweighted", "IPW Weighted")),
    covariate = dplyr::recode(covariate,
      lgdppc_lag = "Lagged Log GDP pc",
      education  = "Education",
      trade      = "Trade Openness"
    )
  )

fig2 <- ggplot2::ggplot(
  balance_long,
  ggplot2::aes(x = SMD, y = covariate, colour = Type, shape = Type)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_vline(xintercept = 0,    linetype = "solid",  colour = "grey50") +
  ggplot2::geom_vline(xintercept =  0.1, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_vline(xintercept = -0.1, linetype = "dashed", colour = "grey50") +
  ggplot2::scale_colour_manual(values = c("Unweighted" = "tomato", "IPW Weighted" = "steelblue")) +
  ggplot2::labs(
    x       = "Standardised Mean Difference",
    y       = NULL,
    colour  = NULL,
    shape   = NULL,
    title   = "Figure 2: Covariate Balance Before and After IPW Weighting",
    caption = "Dashed lines at |SMD| = 0.1 (acceptable balance threshold)."
  ) +
  theme_replication()

ggplot2::ggsave(
  here::here("output", "figure2_balance_love_plot.pdf"),
  fig2, width = 7, height = 4
)
ggplot2::ggsave(
  here::here("output", "figure2_balance_love_plot.png"),
  fig2, width = 7, height = 4, dpi = 150
)
message("Figure 2 saved.")

# ============================================================
# Figure 3: Point Estimates with 95% CIs Across Estimators
# ============================================================

# Build a data frame of estimates and 95% CIs
est_plot_df <- results_table %>%
  dplyr::mutate(
    ci_lo = Estimate - 1.96 * Std.Error,
    ci_hi = Estimate + 1.96 * Std.Error,
    Estimator = factor(
      Estimator,
      levels = c("IPW / MSM (preferred)", "First Differences",
                 "Two-Way FE", "Pooled OLS")
    )
  )

fig3 <- ggplot2::ggplot(
  est_plot_df,
  ggplot2::aes(x = Estimate, y = Estimator, xmin = ci_lo, xmax = ci_hi,
               colour = Estimator == "IPW / MSM (preferred)")
) +
  ggplot2::geom_pointrange(size = 0.7) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::scale_colour_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "grey30"),
    guide  = "none"
  ) +
  ggplot2::labs(
    x       = "Estimated Effect of Democracy on Log GDP per Capita",
    y       = NULL,
    title   = "Figure 3: Democracy Effect Estimates Across Estimators",
    caption = "Point estimates with 95% cluster-robust confidence intervals (clustered by country).\nBlue: preferred IPW/MSM estimator."
  ) +
  theme_replication()

ggplot2::ggsave(
  here::here("output", "figure3_estimates_comparison.pdf"),
  fig3, width = 7, height = 4
)
ggplot2::ggsave(
  here::here("output", "figure3_estimates_comparison.png"),
  fig3, width = 7, height = 4, dpi = 150
)
message("Figure 3 saved.")

# ============================================================
# Figure 4: Propensity Score Overlap
# ============================================================

fig4 <- ggplot2::ggplot(
  panel,
  ggplot2::aes(x = ps_denom, fill = factor(democracy), colour = factor(democracy))
) +
  ggplot2::geom_density(alpha = 0.3) +
  ggplot2::scale_fill_manual(
    values = c("0" = "tomato", "1" = "steelblue"),
    labels = c("0" = "Non-democracy (D=0)", "1" = "Democracy (D=1)")
  ) +
  ggplot2::scale_colour_manual(
    values = c("0" = "tomato", "1" = "steelblue"),
    labels = c("0" = "Non-democracy (D=0)", "1" = "Democracy (D=1)")
  ) +
  ggplot2::labs(
    x       = "Estimated Propensity Score (Pr(Democracy = 1 | covariates))",
    y       = "Density",
    fill    = NULL,
    colour  = NULL,
    title   = "Figure 4: Propensity Score Overlap",
    caption = "Good overlap (common support) is needed for IPW to work well."
  ) +
  theme_replication()

ggplot2::ggsave(
  here::here("output", "figure4_ps_overlap.pdf"),
  fig4, width = 7, height = 4
)
ggplot2::ggsave(
  here::here("output", "figure4_ps_overlap.png"),
  fig4, width = 7, height = 4, dpi = 150
)
message("Figure 4 saved.")

message("\nAll figures written to output/")
