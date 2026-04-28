# CI-with-TSCS-Replication

Replication materials for:

> Blackwell, Matthew, and Adam Glynn. 2018. "How to Make Causal Inferences with Time-Series Cross-Sectional Data under Selection on Observables." *American Political Science Review* 112(4): 1067–1082. https://doi.org/10.1017/S0003055418000357

---

## Overview

This repository replicates the analysis in Blackwell & Glynn (2018), which develops methods for causal inference with time-series cross-sectional (TSCS) panel data when treatment assignment depends on observed covariates ("selection on observables"). The paper introduces a framework based on **marginal structural models (MSMs)** and **inverse probability weighting (IPW)** adapted to the TSCS setting, and applies it to a substantive empirical example.

### Key Contributions of the Paper

1. **Identification**: Formalises what "selection on observables" means in a dynamic TSCS context, where current treatment can depend on the history of covariates and past outcomes.
2. **Estimation via IPW**: Constructs stabilised inverse-probability weights from a sequence of propensity-score models (one per period or pooled), then fits weighted regressions.
3. **Marginal Structural Models**: Shows how weighted fixed-effects (or first-difference) regressions recover the average causal effect of a treatment history.
4. **Comparison of estimators**: Demonstrates that naive OLS and standard fixed-effects estimators are biased when there is time-varying confounding, while the IPW estimator is consistent.
5. **Empirical Application**: Applies the method to panel data on democracy and economic growth (following Acemoglu et al. 2008/2019).

---

## Repository Structure

```
CI-with-TSCS-Replication/
├── README.md                  # This file
├── replication.Rmd            # Unified R Markdown document (runs full replication)
├── R/
│   ├── 00_setup.R             # Install and load required packages
│   ├── 01_data_prep.R         # Download / load and clean the data
│   ├── 02_analysis.R          # Core estimation functions (OLS, FE, IPW/MSM)
│   ├── 03_tables.R            # Reproduce paper tables
│   └── 04_figures.R           # Reproduce paper figures
├── data/
│   └── README.md              # Instructions for obtaining the replication data
├── output/                    # Tables (.tex / .csv) and figures (.pdf / .png) written here
└── .github/
    └── workflows/
        └── ci.yml             # GitHub Actions CI: runs replication on every push
```

---

## Getting Started

### 1. Obtain the Data

The replication data are archived on Harvard Dataverse:

> https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QYEOLJ

Download the dataset files and place them in the `data/` directory. See `data/README.md` for details.

The scripts also contain a **built-in simulation fallback**: if the real data files are not present the code generates synthetic TSCS data with the same structure so that all methods and output formats can be exercised and the CI pipeline always passes.

### 2. Install R Dependencies

```r
source("R/00_setup.R")
```

Required packages (all on CRAN):

| Package | Purpose |
|---------|---------|
| `tidyverse` | Data wrangling and `ggplot2` graphics |
| `plm` | Panel linear models (FE, FD, pooled OLS) |
| `ipw` | Inverse probability weighting utilities |
| `sandwich` | HC/HAC robust and cluster-robust standard errors |
| `lmtest` | Coefficient tests with robust SEs |
| `broom` | Tidy model summaries |
| `knitr` / `kableExtra` | Table rendering inside R Markdown |
| `haven` | Read Stata `.dta` replication files |
| `here` | Portable file paths |

### 3. Run the Full Replication

**Option A — Render the R Markdown document (recommended)**

```r
rmarkdown::render("replication.Rmd")
```

This produces `replication.html` (or PDF if LaTeX is installed) containing all tables and figures in order.

**Option B — Run scripts sequentially**

```r
source("R/00_setup.R")
source("R/01_data_prep.R")
source("R/02_analysis.R")
source("R/03_tables.R")
source("R/04_figures.R")
```

Outputs are saved to the `output/` directory.

---

## Continuous Integration

A GitHub Actions workflow (`.github/workflows/ci.yml`) runs the full replication pipeline on every push and pull request. It:

1. Installs R and all package dependencies (with caching via `r-lib/actions`).
2. Runs `Rscript R/00_setup.R` through `Rscript R/04_figures.R`.
3. Uploads the contents of `output/` as a workflow artifact so you can inspect tables and figures without running the code locally.

The simulation fallback ensures the CI always has data to work with.

---

## Methods Implemented

### Notation (following the paper)

- $i = 1, \dots, N$: cross-sectional units (e.g., countries)
- $t = 1, \dots, T$: time periods
- $D_{it} \in \{0,1\}$: binary treatment (e.g., democracy indicator)
- $Y_{it}$: outcome (e.g., log GDP per capita)
- $\mathbf{X}_{it}$: time-varying confounders (lagged outcome, education, trade openness, …)

### 1. Propensity Score / IPW Weights

A pooled logistic regression (or period-specific logit) is estimated:

$$\Pr(D_{it} = 1 \mid \bar{D}_{i,t-1}, \bar{X}_{it}) = \text{logit}^{-1}(\alpha + \boldsymbol{\beta}^\top \mathbf{X}_{it} + \gamma D_{i,t-1})$$

**Stabilised weights** (Robins, Hernán & Brumback 2000):

$$w_{it} = \prod_{s=1}^{t} \frac{\Pr(D_{is} \mid \bar{D}_{i,s-1})}{\Pr(D_{is} \mid \bar{D}_{i,s-1}, \bar{X}_{is})}$$

where the numerator uses a model without time-varying covariates.

### 2. Estimators Compared

| Estimator | Implementation | Bias under time-varying confounding? |
|-----------|---------------|--------------------------------------|
| Pooled OLS | `lm(Y ~ D + X)` | Yes |
| Fixed Effects (within) | `plm(..., model = "within")` | Yes (if covariates are mediators) |
| First Differences | `plm(..., model = "fd")` | Yes (same issue) |
| **IPW / MSM** | Weighted FE with stabilised weights | No (consistent under SOO) |

### 3. Marginal Structural Model

Fit the MSM using weighted least squares on the within-transformed data:

$$Y_{it} - \bar{Y}_i = \delta (D_{it} - \bar{D}_i) + \varepsilon_{it}, \quad \text{weights } w_{it}$$

Standard errors are cluster-robust at the unit level.

---

## Citation

If you use these replication materials, please cite the original paper:

```bibtex
@article{blackwell2018how,
  title   = {How to Make Causal Inferences with Time-Series Cross-Sectional Data
             under Selection on Observables},
  author  = {Blackwell, Matthew and Glynn, Adam N.},
  journal = {American Political Science Review},
  volume  = {112},
  number  = {4},
  pages   = {1067--1082},
  year    = {2018},
  doi     = {10.1017/S0003055418000357}
}
```

---

## License

Replication code is released under the [MIT License](LICENSE). The underlying data follow the terms of the Harvard Dataverse deposit.
