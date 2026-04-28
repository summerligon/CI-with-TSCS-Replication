# ============================================================
# 01_data_prep.R
# Load and clean the Blackwell & Glynn (2018) replication data.
#
# Priority order:
#   1. Real Stata .dta file from Harvard Dataverse in data/
#   2. Simulated synthetic TSCS data (fallback for CI / demo)
# ============================================================

source(here::here("R", "00_setup.R"))

# ------------------------------------------------------------------
# Path to the real replication data (update if filename differs)
# ------------------------------------------------------------------
DATA_FILE <- here::here("data", "panel_data.dta")

# ------------------------------------------------------------------ #
#  Helper: Simulate TSCS data                                         #
#  Mimics the structure of the empirical application:                 #
#    N countries, T periods                                           #
#    Treatment: binary democracy indicator (D_it)                     #
#    Outcome:   log GDP per capita (Y_it)                             #
#    Confounders: education, trade openness                           #
# ------------------------------------------------------------------ #
simulate_tscs_data <- function(N = 100, T_periods = 20, seed = 42) {
  set.seed(seed)
  message("Simulating TSCS data (N=", N, ", T=", T_periods, ") ...")

  # Country-level fixed effects
  alpha_i <- rnorm(N, mean = 8, sd = 1)   # baseline log GDP
  names(alpha_i) <- paste0("C", sprintf("%03d", seq_len(N)))

  # Time fixed effects
  tau_t <- seq(0, 0.3, length.out = T_periods)   # mild time trend

  df_list <- vector("list", N)

  for (i in seq_len(N)) {
    cname  <- names(alpha_i)[i]
    y_lag  <- alpha_i[i]
    d_lag  <- rbinom(1, 1, 0.3)
    rows   <- vector("list", T_periods)

    for (t in seq_len(T_periods)) {
      # Time-varying covariates (confounders)
      educ  <- 40 + 0.5 * t + rnorm(1, sd = 5)   # education (% enrollment)
      trade <- 50 + rnorm(1, sd = 10)             # trade openness

      # Propensity of treatment depends on history
      lp_d <- -1.5 + 0.8 * d_lag + 0.02 * educ + 0.01 * trade + 0.1 * (y_lag - 8)
      prob_d <- plogis(lp_d)
      d_it <- rbinom(1, 1, prob_d)

      # Potential outcomes — true ATT of democracy is +0.3
      y_it <- alpha_i[i] + tau_t[t] +
        0.30 * d_it +        # true causal effect
        0.05 * educ +
        0.02 * trade +
        0.30 * (y_lag - 8) + # AR(1) in outcome
        rnorm(1, sd = 0.3)

      rows[[t]] <- data.frame(
        country      = cname,
        year         = 1990 + t - 1,
        democracy    = d_it,
        lgdppc       = y_it,
        education    = educ,
        trade        = trade,
        lgdppc_lag   = y_lag,
        democracy_lag = d_lag
      )
      y_lag <- y_it
      d_lag <- d_it
    }
    df_list[[i]] <- dplyr::bind_rows(rows)
  }

  panel <- dplyr::bind_rows(df_list)
  message("Simulated ", nrow(panel), " observations (", N, " units x ", T_periods, " periods).")
  panel
}

# ------------------------------------------------------------------ #
#  Load data                                                          #
# ------------------------------------------------------------------ #
if (file.exists(DATA_FILE)) {
  message("Loading real replication data from: ", DATA_FILE)
  raw <- haven::read_dta(DATA_FILE)

  # --- Harmonise column names to what the rest of the scripts expect ---
  # Adjust the variable_map below to match the actual variable names in your .dta file.
  # The left side is the standard name used in all downstream scripts;
  # the right side is a vector of possible names in the replication data file —
  # the first match found is used.
  variable_map <- list(
    country   = c("country", "ccode", "countryname", "wbcode"),
    year      = c("year", "yr"),
    democracy = c("dem", "democracy", "polity2_bin", "fhfree"),
    lgdppc    = c("lgdppc", "lrgdpch", "log_gdppc", "lgdp"),
    education = c("education", "educ", "school", "sec_school"),
    trade     = c("trade", "tradewb", "openness")
  )

  # Build a rename vector: new_name = old_name (only for cols that need renaming)
  raw_names   <- names(raw)
  rename_vec  <- vapply(names(variable_map), function(std) {
    candidates <- variable_map[[std]]
    match      <- intersect(candidates, raw_names)
    if (length(match) == 0 || match[1] == std) return(NA_character_)
    match[1]                      # first candidate present in data
  }, character(1))
  rename_vec <- rename_vec[!is.na(rename_vec)]  # drop already-correct names

  panel_raw <- if (length(rename_vec) > 0) {
    dplyr::rename(raw, !!!stats::setNames(rename_vec, names(rename_vec)))
  } else {
    raw
  }
  USE_SIMULATED <- FALSE
} else {
  message("Real data file not found at: ", DATA_FILE)
  message("Falling back to simulated data. To use the real data, see data/README.md.")
  panel_raw <- simulate_tscs_data(N = 100, T_periods = 20)
  USE_SIMULATED <- TRUE
}

# ------------------------------------------------------------------ #
#  Clean and construct analysis variables                             #
# ------------------------------------------------------------------ #
panel <- panel_raw %>%
  dplyr::arrange(country, year) %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(
    # Ensure lags exist (use pre-computed if available, else create)
    lgdppc_lag    = dplyr::lag(lgdppc,    1),
    democracy_lag = dplyr::lag(democracy, 1)
  ) %>%
  dplyr::ungroup() %>%
  # Drop first period per country (missing lag)
  dplyr::filter(!is.na(lgdppc_lag), !is.na(democracy_lag)) %>%
  # Ensure treatment is integer 0/1
  dplyr::mutate(democracy = as.integer(as.numeric(democracy)))

message("Final analytic sample: ", nrow(panel), " observations, ",
        dplyr::n_distinct(panel$country), " countries, ",
        dplyr::n_distinct(panel$year), " years.")

# Quick summary
print(summary(panel[, c("lgdppc", "democracy", "education", "trade")]))
