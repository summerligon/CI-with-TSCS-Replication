# advertising-sim.R
# ---------------------------------------------------------------
# Extension of Blackwell & Glynn (2018) to a marketing application
# Simulates TSCS data for N products x T time periods
# Treatment: continuous ad spending (log)
# Outcome: sales next period
# Time-varying confounder: price
# Two conditions:
#   - Endogenous: price affected by past ad spending (treatment-induced confounding)
#   - Exogenous:  price NOT affected by past ad spending (no confounding)
# ---------------------------------------------------------------

library(mvtnorm)
library(mgcv)
library(sandwich)
library(lmtest)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(plm)
library(foreign)
library(reshape)

# Auto-detect repo
possible_repos <- c(
  "/workspaces/Replication-of-Glynn-CI-TSCS",
  "/workspaces/CI-with-TSCS-Replication"
)
repo <- possible_repos[dir.exists(possible_repos)][1]
message("Using repo: ", repo)

panel_utils_paths <- c(
  file.path(repo, "R", "panel-utils.R"),
  file.path(repo, "code", "panel-utils.R"),
  file.path(repo, "R", "code", "code", "panel-utils.R")
)
panel_utils <- panel_utils_paths[file.exists(panel_utils_paths)][1]
if (is.na(panel_utils)) stop("Cannot find panel-utils.R!")
source(panel_utils)

dir.create(file.path(repo, "output"), showWarnings = FALSE)
set.seed(20240429)

# True effect parameters
true.beta <- c(0.40, 0.25, 0.10, 0.05)

# ---------------------------------------------------------------
# 1. DATA GENERATING PROCESS
# endogenous = TRUE:  price affected by past spending (confounding)
# endogenous = FALSE: price NOT affected by past spending (no confounding)
# confound.strength: how strongly past spending affects price
# ---------------------------------------------------------------
sim.adv.data <- function(N = 20, TT = 50, seed = NULL,
                         endogenous = TRUE,
                         confound.strength = -0.7) {
  if (!is.null(seed)) set.seed(seed)
  prod.fe <- rnorm(N, 0, 0.5)

  df <- data.frame(
    product = rep(1:N, each = TT),
    time    = rep(1:TT, times = N),
    sales   = NA, adspend = NA, price = NA, season = NA
  )

  for (i in 1:N) {
    idx     <- which(df$product == i)
    adspend <- rep(NA, TT)
    price   <- rep(NA, TT)
    sales   <- rep(NA, TT)
    season  <- sin(2 * pi * (1:TT) / 12)

    adspend[1:3] <- rnorm(3, 2 + prod.fe[i], 0.3)
    price[1:3]   <- rnorm(3, 5 - 0.2 * adspend[1:3], 0.3)
    sales[1:3]   <- rnorm(3, 3 + 0.3 * adspend[1:3], 0.5)

    for (t in 4:TT) {
      # Price: endogenous = affected by past spending, exogenous = not
      price.effect <- if (endogenous) confound.strength * adspend[t-1] else 0
      price[t] <- 4.0 +
        price.effect +
        0.4 * price[t-1] +
        prod.fe[i] * 0.5 +
        rnorm(1, 0, 0.3)

      # Ad spending: affected by past spending, price, seasonality
      adspend[t] <- 1.5 +
        0.4 * adspend[t-1] +        # reduced from 0.5 to reduce collinearity
        0.2 * price[t-1] +
        0.3 * season[t] +
        prod.fe[i] +
        rnorm(1, 0, 0.4)

      # Sales: true causal effects of ad spending lags
      if (t < TT) {
        sales[t+1] <- 2.0 +
          true.beta[1] * adspend[t] +
          true.beta[2] * adspend[t-1] +
          true.beta[3] * adspend[t-2] +
          true.beta[4] * adspend[t-3] +
          (-0.3) * price[t] +
          0.3 * season[t] +
          prod.fe[i] * 0.4 +
          rnorm(1, 0, 0.4)
      }
    }
    df$adspend[idx] <- adspend
    df$price[idx]   <- price
    df$sales[idx]   <- sales
    df$season[idx]  <- season
  }
  df <- df[complete.cases(df), ]
  return(df)
}

# ---------------------------------------------------------------
# 2. RMSE SIMULATION FUNCTION
# ---------------------------------------------------------------
run.sim <- function(N, TT, n.sims = 200, endogenous = TRUE,
                    confound.strength = -0.7) {
  results <- matrix(NA, nrow = n.sims, ncol = 2,
                    dimnames = list(NULL, c("ADL","SNMM")))
  for (s in 1:n.sims) {
    tryCatch({
      d <- sim.adv.data(N=N, TT=TT, seed=s,
                        endogenous=endogenous,
                        confound.strength=confound.strength)
      d$logspend <- log(d$adspend - min(d$adspend) + 1)
      d$logsales <- log(d$sales   - min(d$sales, na.rm=TRUE) + 1)

      # ADL: includes price as control (biased when endogenous)
      adl <- lm(logsales ~ logspend + l(logspend, product) +
                  price + l(price, product) +
                  season + as.factor(product), data=d)

      # SNMM: blip-down one step
      snm0 <- lm(logsales ~ logspend + l(logspend, product) +
                   price + l(price, product) +
                   season + as.factor(product), data=d)
      snm1 <- lm(I(logsales - coef(snm0)["logspend"] * logspend) ~
                   l(logspend, product) + l(price, product) +
                   season + as.factor(product), data=d)

      results[s,] <- c(coef(adl)["l(logspend, product)"],
                       coef(snm1)["l(logspend, product)"])
    }, error = function(e) NULL)
  }
  apply(results, 2, function(x) sqrt(mean((x - true.beta[2])^2, na.rm=TRUE)))
}

# ---------------------------------------------------------------
# 3. RUN SIMULATIONS: endogenous vs exogenous x N x T
# ---------------------------------------------------------------
message("\nRunning RMSE simulations...")
N.vals  <- c(10, 20, 50, 100)
TT.vals <- c(20, 50)
n.sims  <- 200

sim.results <- expand.grid(N=N.vals, TT=TT.vals,
                           endogenous=c(TRUE, FALSE))
sim.results$ADL  <- NA
sim.results$SNMM <- NA

for (i in 1:nrow(sim.results)) {
  N   <- sim.results$N[i]
  TT  <- sim.results$TT[i]
  endo <- sim.results$endogenous[i]
  message("  N=", N, ", T=", TT,
          ", ", ifelse(endo, "Endogenous", "Exogenous"), "...")
  rmse <- run.sim(N=N, TT=TT, n.sims=n.sims, endogenous=endo)
  sim.results$ADL[i]  <- rmse["ADL"]
  sim.results$SNMM[i] <- rmse["SNMM"]
}

message("Simulation complete!")
print(sim.results)

# ---------------------------------------------------------------
# 4. FIGURE: 2x2 RMSE plot mirroring Blackwell Figure 5
# Left column:  Price endogenous (treatment-induced confounding)
# Right column: Price exogenous  (no confounding)
# Top row:    T = 20
# Bottom row: T = 50
# ---------------------------------------------------------------
message("\nGenerating Figure 5-style RMSE plot...")

cols <- c("ADL"  = "#E41A1C",
          "SNMM" = "#377EB8")
pchs <- c("ADL"  = 17,
          "SNMM" = 19)
ltys <- c("ADL"  = 2,
          "SNMM" = 1)

pdf(file.path(repo, "output", "fig-adv-rmse.pdf"),
    height = 8, width = 10)
par(mfrow = c(2,2), mar = c(4,4,3,1))

for (tt in TT.vals) {
  for (endo in c(TRUE, FALSE)) {
    sub <- sim.results[sim.results$TT == tt &
                       sim.results$endogenous == endo, ]
    sub <- sub[order(sub$N), ]

    ymax <- max(sim.results[,c("ADL","SNMM")]) * 1.1

    title <- paste0(
      ifelse(endo, "Price endogenous", "Price exogenous"),
      " (T = ", tt, ")"
    )

    plot(NULL, xlim=c(8,105), ylim=c(0, ymax),
         xlab="Sample Size (N)", ylab="RMSE",
         main=title, bty="n", las=1, xaxt="n")
    axis(1, at=N.vals)
    abline(h=0, col="grey80")

    for (est in c("ADL","SNMM")) {
      lines(sub$N, sub[[est]],
            col=cols[est], lty=ltys[est])
      points(sub$N, sub[[est]],
             col=cols[est], pch=pchs[est], cex=1.2)
    }

    legend("topright",
           legend = c("ADL", "SNMM"),
           col    = cols,
           pch    = pchs,
           lty    = ltys,
           bty    = "n", cex = 0.9)
  }
}
dev.off()
message("Saved: fig-adv-rmse.pdf")

# ---------------------------------------------------------------
# 5. MAIN DATASET + IRF/SRF FIGURE (endogenous condition)
# ---------------------------------------------------------------
message("\nGenerating IRF/SRF on main dataset...")
adv.data <- sim.adv.data(N=20, TT=50, seed=20240429,
                         endogenous=TRUE, confound.strength=-0.7)
adv.data$logspend <- log(adv.data$adspend - min(adv.data$adspend) + 1)
adv.data$logsales <- log(adv.data$sales   - min(adv.data$sales, na.rm=TRUE) + 1)

mod0 <- lm(logsales ~ logspend + l(logspend, product) +
             price + l(price, product) +
             season + as.factor(product), data=adv.data)
lag1 <- lm(I(logsales - coef(mod0)["logspend"]*logspend) ~
             l(logspend, product) + l(price, product) +
             season + as.factor(product), data=adv.data)
lag2 <- lm(I(logsales -
               coef(mod0)["logspend"]*logspend -
               coef(lag1)["l(logspend, product)"]*l(logspend, product)) ~
             l(logspend, product, 2) + l(price, product, 2) +
             season + as.factor(product), data=adv.data)
lag3 <- lm(I(logsales -
               coef(mod0)["logspend"]*logspend -
               coef(lag1)["l(logspend, product)"]*l(logspend, product) -
               coef(lag2)["l(logspend, product, 2)"]*l(logspend, product, 2)) ~
             l(logspend, product, 3) + l(price, product, 3) +
             season + as.factor(product), data=adv.data)

mods    <- list(mod0, lag1, lag2, lag3)
vcv     <- snmm.var(mods=mods, blip.vars=2,
                    data=adv.data, unit.var="product")
ses     <- sqrt(diag(vcv))
ii      <- c(0, cumsum(sapply(mods, function(x) length(coef(x))))[1:3])
eff.est <- sapply(mods, function(x) coef(x)[2])
eff.ses <- ses[2 + ii]
snmm.srf     <- cumsum(eff.est)
snmm.srf.ses <- c(eff.ses[1],
  sqrt(sum(vcv[ii[1:2]+2, ii[1:2]+2])),
  sqrt(sum(vcv[ii[1:3]+2, ii[1:3]+2])),
  sqrt(sum(vcv[ii[1:4]+2, ii[1:4]+2])))

# ADL comparison
robust.se <- function(fm, clvar) {
  x <- eval(fm$call$data, envir=parent.frame())
  cluster <- x[names(predict(fm)), clvar]
  M <- length(unique(cluster)); N <- length(cluster); K <- dim(vcov(fm))[1]
  dfc <- (M/(M-1)) * ((N-1)/(N-K))
  uj  <- apply(estfun(fm), 2, function(x) tapply(x, cluster, sum))
  dfc * sandwich(fm, meat=crossprod(uj)/N)
}
aa <- coef(mod0)["price"]; b1 <- coef(mod0)["logspend"]
b2 <- coef(mod0)["l(logspend, product)"]
nms <- c("price","logspend","l(logspend, product)")
vcovCL  <- robust.se(mod0, "product")
adl.vcv <- vcovCL[nms, nms]
h1 <- c(b1, aa, 1); h2 <- c(2*aa*b1, aa^2, aa); h3 <- c(3*aa^2*b1, aa^3, aa^2)
adl.est <- c(b1, aa*b1+b2, aa^2*b1+aa*b2, aa^3*b1+aa^2*b2)
adl.ses <- c(eff.ses[1],
             sqrt(t(h1)%*%adl.vcv%*%h1),
             sqrt(t(h2)%*%adl.vcv%*%h2),
             sqrt(t(h3)%*%adl.vcv%*%h3))
adl.srf <- cumsum(adl.est)

message("SNMM estimates: ", paste(round(eff.est,3), collapse=", "))
message("ADL estimates:  ", paste(round(adl.est,3), collapse=", "))
message("True effects:   ", paste(true.beta, collapse=", "))

pdf(file.path(repo, "output", "fig-adv-irf.pdf"), height=4, width=7)
par(mfrow=c(1,2))
plot(x=0:3-0.1, y=eff.est, ylim=c(-0.3,0.7), xlim=c(-0.3,3.3),
     pch=19, col="dodgerblue", las=1, xaxt="n", bty="n",
     xlab="Lag", ylab="Effect on Sales",
     main="Impulse Response Function\n(Price Endogenous)")
axis(1, at=0:3)
points(x=0:3+0.1, adl.est, pch=17, col="indianred")
segments(x0=0:3-0.1, y0=eff.est-1.96*eff.ses,
         y1=eff.est+1.96*eff.ses, col="dodgerblue")
segments(x0=0:3+0.1, y0=adl.est-1.96*adl.ses,
         y1=adl.est+1.96*adl.ses, col="indianred")
abline(h=0, col="grey70")
points(x=0:3, y=true.beta, pch=15, col="darkgreen")
legend("topright", legend=c("SNMM","ADL","Truth"),
       col=c("dodgerblue","indianred","darkgreen"),
       pch=c(19,17,15), bty="n", cex=0.8)

plot(x=0:3+0.1, y=adl.srf, ylim=c(-0.3,1.3), xlim=c(-0.3,3.3),
     pch=17, col="indianred", bty="n", xaxt="n",
     xlab="Lag", ylab="Cumulative Effect",
     main="Step Response Function\n(Price Endogenous)")
axis(1, at=0:3)
abline(h=0, col="grey70")
segments(x0=0:3+0.1, y0=adl.srf-1.96*adl.ses,
         y1=adl.srf+1.96*adl.ses, col="indianred")
points(x=0:3-0.1, y=snmm.srf, pch=19, col="dodgerblue")
segments(x0=0:3-0.1, y0=snmm.srf-1.96*snmm.srf.ses,
         y1=snmm.srf+1.96*snmm.srf.ses, col="dodgerblue")
points(x=0:3, y=cumsum(true.beta), pch=15, col="darkgreen")
legend("topright", legend=c("SNMM","ADL","Truth"),
       col=c("dodgerblue","indianred","darkgreen"),
       pch=c(19,17,15), bty="n", cex=0.8)
dev.off()
message("Saved: fig-adv-irf.pdf")

message("\n=== Advertising simulation complete! ===")
message("Outputs in output/:")
message("  fig-adv-irf.pdf  -- IRF and SRF")
message("  fig-adv-rmse.pdf -- RMSE 2x2 (endogenous vs exogenous x T)")
