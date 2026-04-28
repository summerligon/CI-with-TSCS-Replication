# Data Instructions

## Replication Data from Harvard Dataverse

The replication data for Blackwell & Glynn (2018) are archived at:

> **Harvard Dataverse**: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QYEOLJ

### Steps to Download

1. Visit the Dataverse link above.
2. Click **Access Dataset** → **Original Format ZIP** to download all files, or download individual files.
3. Place the downloaded files in this `data/` directory.

The main file you need is typically a Stata `.dta` file (e.g., `panel_data.dta`). The `R/01_data_prep.R` script looks for:

```
data/panel_data.dta
```

If the file is named differently, update the `DATA_FILE` variable at the top of `R/01_data_prep.R`.

---

## Fallback: Simulated Data

If the real data file is **not** present in `data/`, the script `R/01_data_prep.R` automatically generates
a synthetic TSCS dataset that mimics the structure of the paper's data:

- **N = 100** cross-sectional units (countries)
- **T = 20** time periods
- Binary treatment $D_{it}$ (democracy indicator) that depends on lagged covariates
- Continuous outcome $Y_{it}$ (log GDP per capita)
- Time-varying confounders $X1_{it}$, $X2_{it}$ (education, trade openness)

This allows the full analysis pipeline and CI to run without the real data.

---

## Variable Dictionary

| Variable | Description |
|----------|-------------|
| `country` | Country identifier |
| `year` | Year |
| `democracy` | Binary democracy indicator (treatment) |
| `lgdppc` | Log GDP per capita (outcome) |
| `education` | Secondary school enrolment rate |
| `trade` | Trade openness (% GDP) |
| `lgdppc_lag` | Lagged log GDP per capita |
| `democracy_lag` | Lagged democracy indicator |
