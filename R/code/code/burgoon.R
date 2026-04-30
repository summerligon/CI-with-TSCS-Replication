library(here)
dir.create(here::here("output"), showWarnings = FALSE)
library(mvtnorm)
library(foreign)
library(mgcv)
library(plyr)
library(stargazer)
library(geepack)
library(splines)
library(plm)
source(here::here("code/panel-utils.R"))
set.seed(02138)

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

reasonabledata <- read.csv(here::here("data/burgoon.csv"), stringsAsFactors = FALSE)

burg.rep <- lm(sqrt(terrorinclead)~totspendrevlog + govleft + democ  + poplog + govcap + conflict + tradelog  + europe +africa + asia + america + sqrt(terrorinc) + as.factor(year), data=reasonabledata)
summary(burg.rep)
vcovCL <- robust.se(burg.rep, "cow")
coeftest(burg.rep, vcovCL)

## ttotspendrevlog
out <- lm(sqrt(terrorinclead)~totspendrevlog + l(totspendrevlog, cow)  + govleft + democ + l(democ, cow)  + poplog + govcap + conflict + tradelog  + europe +africa + asia + america + sqrt(terrorinc) + sqrt(l(terrorinc, cow)) + as.factor(year), data=reasonabledata)
vcovCL <- robust.se(out, "cow")

mod0.alt <- lm(sqrt(terrorinclead)~totspendrevlog + l(totspendrevlog, cow)  + europe +africa + asia + america + as.factor(year), data=reasonabledata)

lag1 <- lm(I(sqrt(terrorinclead) - coef(out)["totspendrevlog"]*totspendrevlog) ~ l(totspendrevlog, cow)  + l(govleft, cow) + l(democ, cow) + l(poplog, cow) + l(govcap, cow) + l(conflict, cow) + l(tradelog, cow) + europe +africa + asia + america + sqrt(l(terrorinc, cow))+ sqrt(l(terrorinc, cow, 2))+ as.factor(year), data=reasonabledata)


lag2 <- lm(I(sqrt(terrorinclead) - coef(out)["totspendrevlog"]*totspendrevlog - coef(lag1)["l(totspendrevlog, cow)"]*l(totspendrevlog, cow))~l(totspendrevlog, cow, 2)  + l(govleft, cow, 2) + l(democ, cow, 2) + l(democ, cow, 3) + l(poplog, cow, 2) + l(govcap, cow, 2) + l(conflict, cow, 2) + l(tradelog, cow, 2) + europe +africa + asia + america + sqrt(l(terrorinc, cow, 2))+ sqrt(l(terrorinc, cow,3))+ as.factor(year), data=reasonabledata)

lag3 <- lm(I(sqrt(terrorinclead) - coef(out)["totspendrevlog"]*totspendrevlog - coef(lag1)["l(totspendrevlog, cow)"]*l(totspendrevlog, cow) - coef(lag2)["l(totspendrevlog, cow, 2)"]*l(totspendrevlog, cow, 2))~l(totspendrevlog, cow, 3)  + l(govleft, cow, 3) + l(democ, cow, 3) + l(democ, cow, 4) + l(poplog, cow, 3) + l(govcap, cow, 3) + l(conflict, cow, 3) + l(tradelog, cow, 3) + europe +africa + asia + america + sqrt(l(terrorinc, cow, 3))+ sqrt(l(terrorinc, cow, 4))+ as.factor(year), data=reasonabledata)


lag4 <- lm(I(sqrt(terrorinclead) - coef(out)["totspendrevlog"]*totspendrevlog - coef(lag1)["l(totspendrevlog, cow)"]*l(totspendrevlog, cow) - coef(lag2)["l(totspendrevlog, cow, 2)"]*l(totspendrevlog, cow, 2) - coef(lag3)["l(totspendrevlog, cow, 3)"]*l(totspendrevlog, cow, 3))~l(totspendrevlog, cow, 4)  + l(govleft, cow, 4) + l(democ, cow, 4) + l(democ, cow, 5) + l(poplog, cow, 4) + l(govcap, cow, 4) + l(conflict, cow, 4) + l(tradelog, cow, 4) + europe +africa + asia + america + sqrt(l(terrorinc, cow, 4))+ sqrt(l(terrorinc, cow, 5))+ as.factor(year), data=reasonabledata)

lag5 <- lm(I(sqrt(terrorinclead) - coef(out)["totspendrevlog"]*totspendrevlog - coef(lag1)["l(totspendrevlog, cow)"]*l(totspendrevlog, cow) - coef(lag2)["l(totspendrevlog, cow, 2)"]*l(totspendrevlog, cow, 2) - coef(lag3)["l(totspendrevlog, cow, 3)"]*l(totspendrevlog, cow, 3)- coef(lag4)["l(totspendrevlog, cow, 4)"]*l(totspendrevlog, cow, 4))~l(totspendrevlog, cow, 5)  + l(govleft, cow, 5) + l(democ, cow, 5) + l(democ, cow, 6) + l(poplog, cow, 5) + l(govcap, cow, 5) + l(conflict, cow, 5) + l(tradelog, cow, 5) + europe +africa + asia + america + sqrt(l(terrorinc, cow, 5))+sqrt(l(terrorinc, cow, 6))+ as.factor(year), data=reasonabledata)

lag6 <- lm(I(sqrt(terrorinclead) - coef(out)["totspendrevlog"]*totspendrevlog - coef(lag1)["l(totspendrevlog, cow)"]*l(totspendrevlog, cow) - coef(lag2)["l(totspendrevlog, cow, 2)"]*l(totspendrevlog, cow, 2) - coef(lag3)["l(totspendrevlog, cow, 3)"]*l(totspendrevlog, cow, 3)- coef(lag4)["l(totspendrevlog, cow, 4)"]*l(totspendrevlog, cow, 4) - coef(lag5)["l(totspendrevlog, cow, 5)"]*l(totspendrevlog, cow, 5))~l(totspendrevlog, cow, 6)  + l(govleft, cow, 6) + l(democ, cow, 6) + l(democ, cow, 7) + l(poplog, cow, 6) + l(govcap, cow, 6) + l(conflict, cow, 6) + l(tradelog, cow, 6) + europe +africa + asia + america + sqrt(l(terrorinc, cow, 6))+ sqrt(l(terrorinc, cow, 7))+ as.factor(year), data=reasonabledata)


aa <- coef(out)["sqrt(terrorinc)"]
b1 <- coef(out)["totspendrevlog"]
b2 <- coef(out)["l(totspendrevlog, cow)"]
nms <- c("sqrt(terrorinc)", "totspendrevlog", "l(totspendrevlog, cow)")
adl.vcv <- vcovCL[nms, nms]


adl.l1 <- aa * b1 + b2
adl.s1 <- b1 + aa * b1 + b2
h1 <- c(b1, aa, 1)
s.h1 <- c(b1, 1 + aa, 1)
adl.l1.se <- sqrt(t(h1) %*% adl.vcv %*% h1)
adl.s1.se <- sqrt(t(s.h1) %*% adl.vcv %*% s.h1)

adl.l2 <- aa^2 * b1 + aa * b2
adl.s2 <- adl.l2 + adl.s1
h2 <- c(2 * aa * b1, aa^2, aa)
s.h2 <- h2 + s.h1
adl.l2.se <- sqrt(t(h2) %*% adl.vcv %*% h2)
adl.s2.se <- sqrt(t(s.h2) %*% adl.vcv %*% s.h2)

adl.l3 <- aa^3 * b1 + aa^2 * b2
adl.s3 <- adl.l3 + adl.s2
h3 <- c(3 * aa^2 * b1, aa^3, aa^2)
s.h3 <- h3 + s.h2
adl.l3.se <- sqrt(t(h3) %*% adl.vcv %*% h3)
adl.s3.se <- sqrt(t(s.h3) %*% adl.vcv %*% s.h3)

adl.l4 <- aa^4 * b1 + aa^3 * b2
adl.s4 <- adl.l4 + adl.s3
h4 <- c(4 * aa^3 * b1, aa^4, aa^3)
s.h4 <- h4 + s.h3
adl.l4.se <- sqrt(t(h4) %*% adl.vcv %*% h4)
adl.s4.se <- sqrt(t(s.h4) %*% adl.vcv %*% s.h4)


mods <- list(out, lag1, lag2, lag3, lag4)
vcv <- snmm.var(mods = mods, blip.vars = 2, data = reasonabledata, unit.var = "cow")
ses <- sqrt(diag(vcv))

ii <- c(0,cumsum(lapply(mods, function(x) length(coef(x))))[1:(length(mods)-1)])
eff.est <- unlist(lapply(mods, function(x) coef(x)[2]))
eff.ses <- ses[2+ii]

## SRF for ADL model
adl.est <- c(b1,adl.l1,adl.l2,adl.l3,adl.l4)
adl.ses <- c(eff.ses[1],adl.l1.se,adl.l2.se,adl.l3.se,adl.l4.se)
adl.srf <- c(b1, adl.s1, adl.s2, adl.s3, adl.s4)
adl.srf.ses <- c(eff.ses[1], adl.s1.se, adl.s2.se, adl.s3.se, adl.s4.se)

## SRF for SNMM
snmm.srf <- cumsum(eff.est)
snmm.srf.ses <- c(eff.ses[1],
                  sqrt(sum(vcv[ii[1:2] + 2, ii[1:2] +2])),
                  sqrt(sum(vcv[ii[1:3] + 2, ii[1:3] +2])),
                  sqrt(sum(vcv[ii[1:4] + 2, ii[1:4] +2])),
                  sqrt(sum(vcv[ii[1:5] + 2, ii[1:5] +2])))

snmm.srf/adl.srf

cairo_pdf(here::here("output/fig6-burgoon.pdf"), family = "Minion Pro", height = 4, width = 7, pointsize = 11)
par(mfrow = c(1,2))
plot(x= 0:4 - 0.15, y=eff.est, ylim = c(-0.75, .5), xlim = c(-0.5, 4.25), pch = 19, col = "dodgerblue", las = 1, xaxt = "n", bty = "n", xlab = "Lag of Effect", ylab = "Effect of Govt Spending on Terrorist Incidents", main = "Impulse Response Function")
#legend(x="topright", legend = c("SNMM", "ADL"), col = c("dodgerblue", "indianred"), pch = c(19, 17), bty = "n", lwd = 1)
axis(1, at = 0:4)
points(x=0:4 + 0.15, adl.est, pch = 17, col = "indianred")
segments(x0=0:4 + 0.15, y0 = adl.est - 1.96 * adl.ses, y1 = adl.est + 1.96 * adl.ses, col = "indianred")
text(1 - 0.15, eff.est[2], "SNMM", pos = 2, col = "dodgerblue")
text(1 + 0.15, adl.est[2], "ADL", pos = 4, col = "indianred")
#points(x=0:4 + 0.2, coef(mod0.alt)[2:6], pch = 15, col = "orange")
#segments(x0=0:4 + 0.2, y0=confint(mod0.alt, level = 0.95)[2:6,1], y1=confint(mod0.alt, level = 0.95)[2:6,2], col = "orange")
segments(x0=0:4 - 0.15, y0=eff.est - 1.96*eff.ses, y1 = eff.est + 1.96 * eff.ses, col = "dodgerblue")
abline(h=0, col = "grey70")
plot(x = 0:4 + 0.15, y = adl.srf, ylim = c(-0.75, 0.5), xlim = c(-0.5, 4.25), pch = 17, col = "indianred", xlab = "Lag of Effect", ylab = "Cumulative Effect of Gov't Spending", bty = "n", main = "Step Response Function", las = 1)
abline(h = 0, col = "grey70")
#legend(x="topright", legend = c("SNMM", "ADL"), col = c("dodgerblue", "indianred"), pch = c(19, 17), lty = c(1,2), bty = "n", lwd = 1)
text(1 - 0.15, snmm.srf[2], "SNMM", pos = 2, col = "dodgerblue")
text(1 + 0.15, adl.srf[2] + 0.05, "ADL", pos = 4, col = "indianred")
#lines(x = 0:4 + 0.15, y = adl.srf, lty = 2, col = "indianred")
segments(x0=0:4 + 0.15, y0 = adl.srf - 1.96 * adl.srf.ses, y1 = adl.srf + 1.96 * adl.srf.ses, col = "indianred")
points(x = 0:4 -0.15, y = snmm.srf, pch = 19, col = "dodgerblue")
#lines(x = 0:4 -0.15, y = snmm.srf, col = "dodgerblue")
segments(x0=0:4 -0.15, y0=snmm.srf - 1.96*snmm.srf.ses, y1 = snmm.srf + 1.96 * snmm.srf.ses, col = "dodgerblue")
dev.off()



## leftist government regression

govleft.reg <- lm(govleft ~ l(totspendrevlog, cow)  + l(govleft, cow) + l(democ, cow) + l(poplog, cow) + l(govcap, cow) + l(conflict, cow) + l(tradelog, cow) + europe +africa + asia + america + sqrt(l(terrorinc, cow))+ sqrt(l(terrorinc, cow, 2))+ as.factor(year), data=reasonabledata)
summary(govleft.reg)



## IPTW

reasonabledata$hi.spend <- 1 * (reasonabledata$totspendrevper > 25)
reasonabledata$hi.spend.l1 <- pan.lag(reasonabledata$hi.spend, reasonabledata$cow)
reasonabledata$hi.spend.l2 <- pan.lag(reasonabledata$hi.spend.l1, reasonabledata$cow)
reasonabledata$hi.spend.l3 <- pan.lag(reasonabledata$hi.spend.l2, reasonabledata$cow)
reasonabledata$hi.spend.sum <- pan.sum(reasonabledata$hi.spend, reasonabledata$cow)
reasonabledata$hi.spend.sum.l1 <- pan.sum(reasonabledata$hi.spend.l1, reasonabledata$cow)
reasonabledata$hi.spend.sum.l2 <- pan.sum(reasonabledata$hi.spend.l2, reasonabledata$cow)
reasonabledata$hi.spend.sum.l3 <- pan.sum(reasonabledata$hi.spend.l3, reasonabledata$cow)

burg.ptb <- lm(sqrt(terrorinclead)~ hi.spend + hi.spend.sum.l1  + govleft + democ  + poplog + govcap + conflict + tradelog  + europe + africa + asia + america + sqrt(terrorinc), data=reasonabledata)
coeftest(burg.ptb, robust.se(burg.ptb, "cow"))

burg.omit <- lm(sqrt(terrorinclead)~ hi.spend + hi.spend.sum.l1  + europe + africa + asia + america + sqrt(terrorinc), data=reasonabledata)
coeftest(burg.omit, robust.se(burg.omit, "cow"))

spend.psmod <- as.formula("hi.spend ~ hi.spend.l1 + hi.spend.l2 + hi.spend.sum.l3 + govleft + democ  + poplog + govcap + conflict + tradelog  + europe + africa + asia + america + sqrt(terrorinc)")
spend.marmod <- as.formula("hi.spend ~ hi.spend.l1 +  hi.spend.l2  + hi.spend.sum.l3 + europe + africa + asia + america")

t.wmod <- iptw(denominator = spend.psmod, numerator = spend.marmod,
                    data = reasonabledata, id = "cow", time = "year",
                    family = "binomial")
reasonabledata$sw <- t.wmod$sw
reasonabledata$sw[which(reasonabledata$sw > 10)] <- 10

burg.rep.data <- model.frame(terrorinclead ~ hi.spend + hi.spend.sum.l1  + europe + africa + asia + america + year + sw + cow, data = reasonabledata)
burg.iptw <- geeglm(sqrt(terrorinclead) ~ hi.spend + hi.spend.sum.l1 + europe + africa + asia + america, data = burg.rep.data, id = cow, waves = year, corstr = "independence", weights = sw)
summary(burg.iptw)

burg.iptw <- lm(sqrt(terrorinclead) ~ hi.spend + hi.spend.sum.l1 + europe + africa + asia + america, data = reasonabledata, weights = sw)
coeftest(burg.iptw, robust.se(burg.iptw, "cow"))

boots <- 1000
bhold <- matrix(NA, ncol = 3, nrow = boots)
rowlist <- tapply(1:nrow(reasonabledata), reasonabledata$cow, function(x) x)
n.cows <- length(unique(reasonabledata$cow))
for (b in 1:boots) {
  star <- sample(1:n.cows, size = n.cows, replace = TRUE)
  dstar <- reasonabledata[unlist(rowlist[star]),]
  burg.ptb.star <- lm(sqrt(terrorinclead)~ hi.spend + hi.spend.sum.l1  + govleft + democ  + poplog + govcap + conflict + tradelog  + europe + africa + asia + america + sqrt(terrorinc), data=dstar)
  burg.omit.star <- lm(sqrt(terrorinclead)~ hi.spend + hi.spend.sum.l1  + europe + africa + asia + america + sqrt(terrorinc), data=dstar)

  t.wmod <- iptw(denominator = spend.psmod, numerator = spend.marmod,
                    data = dstar, id = "cow", time = "year",
                    family = "binomial")
  dstar$sw <- t.wmod$sw
  dstar$sw[which(dstar$sw > 10)] <- 10

  burg.iptw.star <- lm(sqrt(terrorinclead) ~ hi.spend + hi.spend.sum.l1 + europe + africa + asia + america, data = dstar, weights = sw)

  bhold[b,1] <- coef(burg.ptb.star)["hi.spend.sum.l1"]
  bhold[b,2] <- coef(burg.omit.star)["hi.spend.sum.l1"]
  bhold[b,3] <- coef(burg.iptw.star)["hi.spend.sum.l1"]
}

b.ests <- c(coef(burg.ptb)['hi.spend.sum.l1'], coef(burg.omit)['hi.spend.sum.l1'],
            coef(burg.iptw)['hi.spend.sum.l1'])
ci.hi <- b.ests + 1.96 * apply(bhold, 2, sd)
ci.lo <- b.ests - 1.96 * apply(bhold, 2, sd)

cairo_pdf(filename = here::here("output/fig7-burgoon-iptw.pdf"), family = "Minion Pro", height = 5, width = 5, pointsize = 11)
plot(x = coef(burg.ptb)['hi.spend.sum.l1'], y = 6, ylim = c(2, 6), xlim = c(-.1,0.1), pch = 19, yaxt = "n", bty = "n", ylab = "",
     xlab = "Effect of Cumulative Lagged Welfare Spending on Terrorism")
points(x = coef(burg.omit)['hi.spend.sum.l1'], y = 4, pch = 19)
points(x = coef(burg.iptw)['hi.spend.sum.l1'], y = 2, pch = 19)
abline(v = 0, col = "grey")
segments(x0 = ci.lo[1], x1 = ci.hi[1], y0 = 6)
segments(x0 = ci.lo[2], x1 = ci.hi[2], y0 = 4)
segments(x0 = ci.lo[3], x1 = ci.hi[3], y0 = 2)
text(x = ci.hi[1], y = 6, pos = 4, label = "(a) Control for TVCs")
text(x = ci.hi[2], y = 4, pos = 4, label = "(b) Omit TVCs")
text(x = 0, y = 2, pos = 4, label = "(c) IPTW")
dev.off()
