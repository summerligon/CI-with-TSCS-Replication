# advertising-sim.R
# ---------------------------------------------------------------
# Extension of Blackwell & Glynn (2018) to a marketing application
# Simulates TSCS data for N products x T time periods
# Treatment: continuous ad spending (log)
# Outcome: sales next period
# Time-varying confounder: price (affected by past ad spending)
# Methods: SNMM (blip-down) + Figure 5-style RMSE simulation
# ---------------------------------------------------------------

library(mvtnorm)
library(mgcv)
library(sandwich)
library(lmtest)
library(ggplot2)
library(RColorBrewer)

repo <- "/workspaces/CI-with-TSCS-Replication"
source(file.path(repo, "R", "panel-utils.R"))
dir.create(file.path(repo, "output"), showWarnings = FALSE)
set.seed(20240429)

# ---------------------------------------------------------------
# 1. DATA GENERATING PROCESS
# ---------------------------------------------------------------
# The DGP mirrors Blackwell's setup:
#   - Ad spending (A_t) affects sales (Y_{t+1})
#   - Past ad spending affects price (Z_t) -- treatment-induced confounding
#   - Price (Z_t) also affects current spending and sales
#
# True structural equations:
#   A_t = f(A_{t-1}, Z_{t-1}, product FE, season) + e_A
#   Z_t = g(A_{t-1}, Z_{t-1}, product FE) + e_Z   <- TVC affected by past A
#   Y_{t+1} = h(A_t, A_{t-1}, Z_t, product FE) + e_Y
#
# True lagged effects of ad spending on sales:
#   lag 0: beta0 = 0.4
#   lag 1: beta1 = 0.25
#   lag 2: beta2 = 0.10
#   lag 3: beta3 = 0.05
# ---------------------------------------------------------------

# True effect parameters
true.beta <- c(0.40, 0.25, 0.10, 0.05)  # lag 0,1,2,3

sim.adv.data <- function(N = 20, TT = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Product fixed effects
  prod.fe <- rnorm(N, 0, 0.5)

  # Initialize storage
  n.obs <- N * TT
  df <- data.frame(
    product = rep(1:N, each = TT),
    time    = rep(1:TT, times = N),
    sales   = NA,
    adspend = NA,
    price   = NA,
    season  = NA
  )

  # Simulate panel
  for (i in 1:N) {
    idx <- which(df$product == i)

    # Initialize first 3 periods
    adspend <- rep(NA, TT)
    price   <- rep(NA, TT)
    sales   <- rep(NA, TT)
    season  <- sin(2 * pi * (1:TT) / 12)  # seasonal pattern

    adspend[1:3] <- rnorm(3, 2 + prod.fe[i], 0.3)
    price[1:3]   <- rnorm(3, 5 - 0.2 * adspend[1:3], 0.3)
    sales[1:3]   <- rnorm(3, 3 + 0.3 * adspend[1:3], 0.5)

    for (t in 4:TT) {
      # Price (TVC): affected by PAST ad spending -- treatment-induced confounding
      price[t] <- 4.0 +
        (-0.3) * adspend[t-1] +        # past spending reduces price (promo)
        0.4    * price[t-1] +           # price persistence
        prod.fe[i] * 0.5 +
        rnorm(1, 0, 0.3)

      # Ad spending: affected by past spending, past price, seasonality
      adspend[t] <- 1.5 +
        0.5  * adspend[t-1] +           # spending persistence
        0.2  * price[t-1] +             # high price -> more advertising
        0.3  * season[t] +              # seasonal ad campaigns
        prod.fe[i] +
        rnorm(1, 0, 0.3)

      # Sales: true causal effect of ad spending (current + lags)
      # Note: conditioning on price here creates collider bias in OLS
      if (t < TT) {
        sales[t+1] <- 2.0 +
          true.beta[1] * adspend[t] +     # lag 0 effect
          true.beta[2] * adspend[t-1] +   # lag 1 effect
          true.beta[3] * adspend[t-2] +   # lag 2 effect
          true.beta[4] * adspend[t-3] +   # lag 3 effect
          (-0.2) * price[t] +             # price reduces sales
          0.3 * season[t] +
          prod.fe[i] * 0.4 +
          rnorm(1, 0, 0.5)
      }
    }

    df$adspend[idx] <- adspend
    df$price[idx]   <- price
    df$sales[idx]   <- sales
    df$season[idx]  <- season
  }

  # Remove incomplete observations
  df <- df[complete.cases(df), ]
  return(df)
}

# ---------------------------------------------------------------
# 2. GENERATE MAIN DATASET
# ---------------------------------------------------------------
message("Generating advertising TSCS data (N=20, T=50)...")
adv.data <- sim.adv.data(N = 20, TT = 50, seed = 20240429)
adv.data$logspend <- log(adv.data$adspend - min(adv.data$adspend) + 1)
adv.data$logsales <- log(adv.data$sales - min(adv.data$sales, na.rm=TRUE) + 1)

message("Data dimensions: ", nrow(adv.data), " observations")

# ---------------------------------------------------------------
# 3. SNMM ESTIMATION (Blip-Down Algorithm)
# Mirrors burgoon.R exactly
# ---------------------------------------------------------------
message("\nEstimating SNMM via blip-down algorithm...")

# Contemporaneous model (lag 0)
mod0 <- lm(logsales ~ logspend + l(logspend, product) +
             price + l(price, product) +
             season + as.factor(product),
           data = adv.data)

# Blip-down lag 1
lag1 <- lm(I(logsales - coef(mod0)["logspend"] * logspend) ~
             l(logspend, product) +
             l(price, product) +
             season + as.factor(product),
           data = adv.data)

# Blip-down lag 2
lag2 <- lm(I(logsales -
               coef(mod0)["logspend"] * logspend -
               coef(lag1)["l(logspend, product)"] * l(logspend, product)) ~
             l(logspend, product, 2) +
             l(price, product, 2) +
             season + as.factor(product),
           data = adv.data)

# Blip-down lag 3
lag3 <- lm(I(logsales -
               coef(mod0)["logspend"] * logspend -
               coef(lag1)["l(logspend, product)"] * l(logspend, product) -
               coef(lag2)["l(logspend, product, 2)"] * l(logspend, product, 2)) ~
             l(logspend, product, 3) +
             l(price, product, 3) +
             season + as.factor(product),
           data = adv.data)

# SNMM variance (sandwich estimator)
mods <- list(mod0, lag1, lag2, lag3)
vcv  <- snmm.var(mods = mods, blip.vars = 2,
                 data = adv.data, unit.var = "product")
ses  <- sqrt(diag(vcv))

ii       <- c(0, cumsum(lapply(mods, function(x) length(coef(x))))[1:3])
eff.est  <- unlist(lapply(mods, function(x) coef(x)[2]))
eff.ses  <- ses[2 + ii]

# Cumulative (step response)
snmm.srf     <- cumsum(eff.est)
snmm.srf.ses <- c(
  eff.ses[1],
  sqrt(sum(vcv[ii[1:2]+2, ii[1:2]+2])),
  sqrt(sum(vcv[ii[1:3]+2, ii[1:3]+2])),
  sqrt(sum(vcv[ii[1:4]+2, ii[1:4]+2]))
)

# ADL comparison
aa      <- coef(mod0)["price"]
b1      <- coef(mod0)["logspend"]
b2      <- coef(mod0)["l(logspend, product)"]
nms     <- c("price", "logspend", "l(logspend, product)")

# Clustered SEs for ADL
robust.se <- function(fm, clvar) {
  library(sandwich); library(lmtest)
  x       <- eval(fm$call$data, envir = parent.frame())
  cluster <- x[names(predict(fm)), clvar]
  M <- length(unique(cluster)); N <- length(cluster); K <- dim(vcov(fm))[1]
  dfc <- (M/(M-1)) * ((N-1)/(N-K))
  uj  <- apply(estfun(fm), 2, function(x) tapply(x, cluster, sum))
  dfc * sandwich(fm, meat = crossprod(uj)/N)
}
vcovCL  <- robust.se(mod0, "product")
adl.vcv <- vcovCL[nms, nms]

adl.l1 <- aa * b1 + b2
adl.l2 <- aa^2 * b1 + aa * b2
adl.l3 <- aa^3 * b1 + aa^2 * b2

h1 <- c(b1, aa, 1)
h2 <- c(2*aa*b1, aa^2, aa)
h3 <- c(3*aa^2*b1, aa^3, aa^2)

adl.est     <- c(b1, adl.l1, adl.l2, adl.l3)
adl.ses     <- c(eff.ses[1],
                 sqrt(t(h1) %*% adl.vcv %*% h1),
                 sqrt(t(h2) %*% adl.vcv %*% h2),
                 sqrt(t(h3) %*% adl.vcv %*% h3))
adl.srf     <- cumsum(adl.est)
adl.srf.ses <- adl.ses  # simplified

message("SNMM estimates (lags 0-3): ", paste(round(eff.est, 3), collapse=", "))
message("True effects (lags 0-3):   ", paste(true.beta, collapse=", "))

# ---------------------------------------------------------------
# 4. FIGURE A: IRF + SRF (mirrors fig6-burgoon.pdf)
# ---------------------------------------------------------------
message("\nGenerating IRF/SRF figure...")
pdf(file.path(repo, "output", "fig-adv-irf.pdf"), height=4, width=7)
par(mfrow=c(1,2))

# IRF
plot(x=0:3 - 0.1, y=eff.est,
     ylim=c(-0.2, 0.6), xlim=c(-0.3, 3.3),
     pch=19, col="dodgerblue", las=1, xaxt="n", bty="n",
     xlab="Lag of Effect",
     ylab="Effect of Ad Spending on Sales",
     main="Impulse Response Function")
axis(1, at=0:3)
points(x=0:3 + 0.1, adl.est, pch=17, col="indianred")
segments(x0=0:3 - 0.1, y0=eff.est - 1.96*eff.ses,
         y1=eff.est + 1.96*eff.ses, col="dodgerblue")
segments(x0=0:3 + 0.1, y0=adl.est - 1.96*adl.ses,
         y1=adl.est + 1.96*adl.ses, col="indianred")
abline(h=0, col="grey70")
points(x=0:3, y=true.beta, pch=15, col="darkgreen")
text(1 - 0.1, eff.est[2], "SNMM", pos=2, col="dodgerblue", cex=0.8)
text(1 + 0.1, adl.est[2], "ADL",  pos=4, col="indianred",  cex=0.8)
text(1,       true.beta[2],"Truth",pos=3, col="darkgreen",  cex=0.8)

# SRF
plot(x=0:3 + 0.1, y=adl.srf,
     ylim=c(-0.2, 1.2), xlim=c(-0.3, 3.3),
     pch=17, col="indianred", bty="n", xaxt="n",
     xlab="Lag of Effect",
     ylab="Cumulative Effect of Ad Spending",
     main="Step Response Function")
axis(1, at=0:3)
abline(h=0, col="grey70")
segments(x0=0:3 + 0.1, y0=adl.srf - 1.96*adl.srf.ses,
         y1=adl.srf + 1.96*adl.srf.ses, col="indianred")
points(x=0:3 - 0.1, y=snmm.srf, pch=19, col="dodgerblue")
segments(x0=0:3 - 0.1, y0=snmm.srf - 1.96*snmm.srf.ses,
         y1=snmm.srf + 1.96*snmm.srf.ses, col="dodgerblue")
points(x=0:3, y=cumsum(true.beta), pch=15, col="darkgreen")
text(1 - 0.1, snmm.srf[2],         "SNMM", pos=2, col="dodgerblue", cex=0.8)
text(1 + 0.1, adl.srf[2],          "ADL",  pos=4, col="indianred",  cex=0.8)
text(1,       cumsum(true.beta)[2], "Truth",pos=3, col="darkgreen",  cex=0.8)
dev.off()
message("Saved: fig-adv-irf.pdf")

# ---------------------------------------------------------------
# 5. FIGURE B: RMSE SIMULATION (mirrors Figure 5 in paper)
# ---------------------------------------------------------------
message("\nRunning RMSE simulation (this takes ~10 min)...")

run.sim <- function(N, TT, n.sims=200) {
  results <- matrix(NA, nrow=n.sims, ncol=2,
                    dimnames=list(NULL, c("ADL","SNMM")))
  for (s in 1:n.sims) {
    tryCatch({
      d <- sim.adv.data(N=N, TT=TT, seed=s)
      d$logspend <- log(d$adspend - min(d$adspend) + 1)
      d$logsales <- log(d$sales   - min(d$sales, na.rm=TRUE) + 1)

      # ADL lag-1 estimate
      adl <- lm(logsales ~ logspend + l(logspend, product) +
                  price + l(price, product) +
                  season + as.factor(product), data=d)
      results[s, "ADL"] <- coef(adl)["l(logspend, product)"]

      # SNMM lag-1 estimate (blip-down one step)
      snm0 <- lm(logsales ~ logspend + l(logspend, product) +
                   price + l(price, product) +
                   season + as.factor(product), data=d)
      snm1 <- lm(I(logsales - coef(snm0)["logspend"]*logspend) ~
                   l(logspend, product) + l(price, product) +
                   season + as.factor(product), data=d)
      results[s, "SNMM"] <- coef(snm1)["l(logspend, product)"]
    }, error=function(e) NULL)
  }
  # Compute RMSE against true lag-1 effect
  truth <- true.beta[2]
  apply(results, 2, function(x) {
    x <- x[!is.na(x)]
    sqrt(mean((x - truth)^2))
  })
}

# Parameter grid (smaller than paper for speed)
N.vals  <- c(10, 20, 50, 100)
TT.vals <- c(20, 50)
n.sims  <- 200

sim.results <- expand.grid(N=N.vals, TT=TT.vals)
sim.results$ADL  <- NA
sim.results$SNMM <- NA

for (i in 1:nrow(sim.results)) {
  N  <- sim.results$N[i]
  TT <- sim.results$TT[i]
  message("  Running N=", N, ", T=", TT, "...")
  rmse <- run.sim(N=N, TT=TT, n.sims=n.sims)
  sim.results$ADL[i]  <- rmse["ADL"]
  sim.results$SNMM[i] <- rmse["SNMM"]
}

message("Simulation complete!")
print(sim.results)

# ---------------------------------------------------------------
# 6. PLOT FIGURE 5-STYLE RMSE FIGURE
# ---------------------------------------------------------------
message("\nGenerating RMSE figure...")

# Reshape for plotting
sim.long <- reshape(sim.results,
                    varying   = c("ADL","SNMM"),
                    v.names   = "RMSE",
                    timevar   = "Estimator",
                    times     = c("ADL","SNMM"),
                    direction = "long")

cols <- c("ADL"  = "#E41A1C",
          "SNMM" = "#377EB8")
pchs <- c("ADL"  = 17,
          "SNMM" = 19)

pdf(file.path(repo, "output", "fig-adv-rmse.pdf"), height=5, width=8)
par(mfrow=c(1,2), mar=c(4,4,3,1))

for (tt in TT.vals) {
  sub <- sim.long[sim.long$TT == tt, ]
  plot(NULL, xlim=c(8,105), ylim=c(0, max(sim.long$RMSE)*1.1),
       xlab="Sample Size (N)", ylab="RMSE",
       main=paste0("T = ", tt, " (Price as TVC)"),
       bty="n", las=1, xaxt="n")
  axis(1, at=N.vals)
  abline(h=0, col="grey80")

  for (est in c("ADL","SNMM")) {
    d <- sub[sub$Estimator == est, ]
    d <- d[order(d$N), ]
    lines(d$N, d$RMSE, col=cols[est], lty=1)
    points(d$N, d$RMSE, col=cols[est], pch=pchs[est])
  }

  if (tt == TT.vals[1]) {
    legend("topright", legend=c("ADL","SNMM"),
           col=cols, pch=pchs, lty=1, bty="n", cex=0.9)
  }
}
dev.off()
message("Saved: fig-adv-rmse.pdf")

message("\n=== Advertising simulation complete! ===")
message("Outputs saved to output/:")
message("  fig-adv-irf.pdf  -- IRF and SRF (mirrors Figure 6)")
message("  fig-adv-rmse.pdf -- RMSE simulation (mirrors Figure 5)")
message("\nTrue lag effects: ", paste(true.beta, collapse=", "))
message("SNMM estimates:   ", paste(round(eff.est, 3), collapse=", "))
message("ADL estimates:    ", paste(round(adl.est, 3), collapse=", "))
