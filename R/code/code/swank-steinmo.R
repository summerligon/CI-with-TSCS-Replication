library(here)
library(foreign)
library(plm)
library(geepack)
library(mgcv)
source(here::here("code/panel-utils.R"))
robust.se <- function(fm, clvar){
    # R-codes (www.r-project.org) for computing
    # clustered-standard errors. Mahmood Arai, Jan 26, 2008.
    # The arguments of the function are:
    # fitted model, cluster1 and cluster2
    # You need to install libraries `sandwich' and `lmtest'
  library(sandwich);library(lmtest);
  x <- eval(fm$call$data, envir = parent.frame())
  if ("polr" %in% class(fm)) {
    require(MASS)
    cluster <- x[rownames(predict(fm, type = "probs")), clvar]
  } else {
    cluster <- x[names(predict(fm)), clvar]
  }
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- dim(vcov(fm))[1]
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  return(vcovCL)
}


swank <- foreign::read.dta(here::here("data/swank.dta"))
is.na(swank) <- swank == -999

pout <- plm(labtx ~ l1capcon + l1trade + l1ltunem + psdebt + l1old + l1labtx + l1grow + l1inflat + l1unemp + l1leftc + l1tcdemc,
            data = swank, subset = coid != 5, model = "pooling")
summary(pout)

swank$btrade.l1 <- 1 * (swank$l1trade >= 53.5)
swank$btradesum.l1 <- pan.sum(swank$btrade.l1, swank$coid)
swank$btrade.l2 <- pan.lag(swank$btrade.l1, swank$coid)
swank$btrade.l3 <- pan.lag(swank$btrade.l2, swank$coid)
swank$btradesum.l3 <- pan.sum(swank$btrade.l3, swank$coid)

swank$bcapcon.l1 <- 1 * (swank$l1capcon >= 3.5)
swank$bcapconsum.l1 <- pan.sum(swank$bcapcon.l1, swank$coid)
swank$bcapcon.l2 <- pan.lag(swank$bcapcon.l1, swank$coid)
swank$bcapcon.l3 <- pan.lag(swank$bcapcon.l2, swank$coid)
swank$bcapconsum.l2 <- pan.sum(swank$bcapcon.l2, swank$coid)
swank$bcapconsum.l3 <- pan.sum(swank$bcapcon.l3, swank$coid)

swank$l1tradesum <- pan.sum(swank$l1trade, swank$coid)
swank$l2trade <- pan.lag(swank$l1trade, swank$coid)
swank$l3trade <- pan.lag(swank$l2trade, swank$coid)
swank$l4trade <- pan.lag(swank$l3trade, swank$coid)
swank$l5trade <- pan.lag(swank$l4trade, swank$coid)
swank$l2tradesum <- pan.sum(swank$l2trade, swank$coid)
swank$l3tradesum <- pan.sum(swank$l3trade, swank$coid)
swank$l4tradesum <- pan.sum(swank$l4trade, swank$coid)
swank$l5tradesum <- pan.sum(swank$l5trade, swank$coid)

pout.bt <- plm(labtx ~ l1capcon + btrade.l1 + l1ltunem + psdebt + l1old + l1labtx + l1grow + l1inflat + l1unemp + l1leftc + l1tcdemc,
            data = swank, subset = coid != 5, model = "pooling")
summary(pout.bt)


swank.tvc.data <- model.frame(labtx ~ l1capcon + btradesum.l1 + l1ltunem + psdebt + l1old + l1labtx + l1grow + l1inflat + l1unemp + l1leftc + l1tcdemc + year + coid, data = swank)
pout.bt.lag <- geeglm(labtx ~ l1capcon + btradesum.l1 + l1ltunem + psdebt + l1old + l1labtx + l1grow + l1inflat + l1unemp + l1leftc + l1tcdemc, data = swank.tvc.data, id = coid, waves = year, corstr = "ar1", subset = coid != 5)
summary(pout.bt.lag)$coefficients

## Weighting model for binary democracy as treatment
trade.psmodel <- as.formula("btrade.l1 ~ btrade.l2 + btradesum.l3 + year + l1capcon+ l1ltunem + psdebt + l1old + l1labtx + l1grow + l1inflat + l1unemp + l1leftc + l1tcdemc")
trade.marmodel <- as.formula("btrade.l1 ~ btrade.l2 +btradesum.l3 + year ")

## This function takes in a formula for the demoninator model and
## numerator model in eq (9) of Blackwell/Glynn and calculates the
## weights in (9). See panel-utils.R for more details.
t.wmod <- iptw(denominator = trade.psmodel, numerator = trade.marmodel,
                    data = swank, id = "coid", time = "year",
                    family = "binomial")
swank$trade.sw <- t.wmod$sw


swank.rep.data <- model.frame(labtx ~ year + btradesum.l1 + trade.sw + coid, data = swank)
swank.rep <- geeglm(labtx ~ year  + btradesum.l1, data = swank.rep.data, id = coid, waves = year, corstr = "ar1", weights = trade.sw)
summary(swank.rep)

swank.rep.data <- model.frame(labtx ~ year + btradesum.l1+ trade.sw + coid, data = swank)
swank.unw <- geeglm(labtx ~ year + btradesum.l1, data = swank.rep.data, id = coid, waves = year, corstr = "ar1")
summary(swank.unw)$coefficients


## cairo_pdf(filename = "si-swank.pdf", family = "Minion Pro", height = 5, width = 5, pointsize = 11)
plot(x = coef(pout.bt.lag)['btradesum.l1'], y = 6, ylim = c(2, 6), xlim = c(-2.1,2.1), pch = 19, yaxt = "n", bty = "n", ylab = "",
     xlab = "Effect of Cumulative Trade Openness on Effective Labor Tax Rate")
points(x = coef(swank.unw)['btradesum.l1'], y = 4, pch = 19)
points(x = coef(swank.rep)['btradesum.l1'], y = 2, pch = 19)
abline(v = 0, col = "grey")
segments(x0 = confint.geelm(pout.bt.lag)['btradesum.l1', 1], x1 = confint.geelm(pout.bt.lag)['btradesum.l1',2], y0 = 6)
segments(x0 = confint.geelm(swank.unw)['btradesum.l1',1], x1 = confint.geelm(swank.unw)['btradesum.l1',2], y0 = 4)
segments(x0 = confint.geelm(swank.rep)['btradesum.l1',1], x1 = confint.geelm(swank.rep)['btradesum.l1',2], y0 = 2)
text(x = confint.geelm(pout.bt.lag)['btradesum.l1',2], y = 6, pos = 4, label = "(a) Control for TVCs")
text(x = confint.geelm(swank.unw)['btradesum.l1',2], y = 4, pos = 4, label = "(b) Omit TVCs")
text(x = 0, y = 2, pos = 2, label = "(c) IPTW")
#dev.off()
