library(here)
## Design for this simulation comes from almirall et al (2013, stat in
## medicine).
library(plm)
library(sandwich)
library(lmtest)
library(RColorBrewer)
source(here::here("code/panel-utils.R"))
inv.logit <- boot::inv.logit
set.seed(92648)

tscs.sim <- function(s, N, T, gamma, ...) {
  y11 <- y10 <- y01 <- y00 <- matrix(NA, nrow = N, ncol = T)
  z <- x <- y <- matrix(NA, nrow = N, ncol = T)

  u <- rnorm(N, 0, sig)

  eps1 <-  0.8 + 0.9 * u
  y00[,1] <- eps1 + rnorm(N, 0, sig)
  y01[,1] <- eps1 + mu2.01 + rnorm(N, 0, sig)
  y10[,1] <- y00[,1]
  y11[,1] <- y00[,1]

  z[,1] <- 1.7 * u + rnorm(N, 0.4, sig)
  p.x <- inv.logit(cbind(1, z[,1]) %*% alpha1)
  x[,1] <- rbinom(N, size = 1, prob = p.x)
  y[,1] <- y01[,1] * x[,1] + y00[,1] * (1 - x[,1])


  for (t in 2:T) {
    eps1 <-  0.8 + 0.9 * u
    e11 <- rnorm(N, 0, sig)
    e10 <- rnorm(N, 0, sig)
    e01 <- rnorm(N, 0, sig)
    e00 <- rnorm(N, 0, sig)

    y11[,t] <- eps1 + mu1.1  + mu2.11 + e11
    y10[,t] <- eps1 + mu1.1  + e10
    y01[,t] <- eps1 + mu2.01 + e01
    y00[,t] <- eps1 + e00

    z1 <- gamma[1] + gamma[2] + 1.7 * u + rnorm(N, 0, sig)
    z0 <- gamma[1] + 1.7 * u + rnorm(N, 0, sig)
    z[,t] <- x[,t-1] * z1 + (1 - x[,t-1]) * z0

    p.x <- cbind(1, z[,t], y[,t-1]) %*% alpha2 + rnorm(N, 0, 1)
    x[,t] <- 1 * (p.x > 0)##rbinom(N, size = 1, prob = p.x)
    y[,t] <- y00[,t] + x[,t-1] * (y10[,t] - y00[,t]) + x[,t] * (y01[,t] - y00[,t]) + x[,t-1] * x[,t] * ((y11[,t] - y01[,t]) - (y10[,t] - y00[,t]))
  }

  ii <- rep(1:N, each = T)
  tt <- rep(1:T, times = N)
  uu <- rep(u, each = T)
  y.vec <- c(t(y))
  x.vec <- c(t(x))
  z.vec <- c(t(z))
  y10.vec <- c(t(y10))
  y00.vec <- c(t(y00))
  y11.vec <- c(t(y11))
  y01.vec <- c(t(y01))

  z.obs <- exp(0.25 * (z.vec - 0.5)^3)


  pan.dat <- data.frame(ii = ii, tt = tt, y = y.vec, x = x.vec, z = z.vec,
                        y00 = y00.vec, y10 = y10.vec, y01 = y01.vec,
                        y11 = y11.vec, u = uu, z.obs = z.obs)
  pan.dat$ly <- pan.lag(pan.dat$y, pan.dat$ii)
  pan.dat$l2y <- pan.lag(pan.dat$y, pan.dat$ii, 2)
  pan.dat$l3y <- pan.lag(pan.dat$y, pan.dat$ii, 3)
  pan.dat$lx <- pan.lag(pan.dat$x, pan.dat$ii)
  pan.dat$l2x <- pan.lag(pan.dat$x, pan.dat$ii, 2)
  pan.dat$l3x <- pan.lag(pan.dat$x, pan.dat$ii, 3)
  pan.dat$lz <- pan.lag(pan.dat$z, pan.dat$ii)
  pan.dat$l2z <- pan.lag(pan.dat$z, pan.dat$ii, 2)
  pan.dat$l3z <- pan.lag(pan.dat$z, pan.dat$ii, 3)
  pan.dat$lzo <- pan.lag(pan.dat$z.obs, pan.dat$ii)
  pan.dat$l2zo <- pan.lag(pan.dat$z.obs, pan.dat$ii, 2)
  pan.dat$l3zo <- pan.lag(pan.dat$z.obs, pan.dat$ii, 3)

  pols <- lm(y ~ x + lx + ly + l2y + z + lz, data = pan.dat)
  pan.dat$ytil <- pan.dat$y - pan.dat$x * coef(pols)["x"]
  pols2 <- lm(ytil ~ l2y + lx + lz, data = pan.dat)
  ## snmm.vcv <- snmm.var(list(pols, pols2), blip.vars = 1, pan.dat, unit.var = "ii")
  ## snmm.se <- sqrt(diag(snmm.vcv)[10])
  pols.alt <- lm(y ~ x + lx + ly, data = pan.dat)

  pols.obs <- lm(y ~ x + lx + ly + l2y + z.obs + lzo, data = pan.dat)
  pan.dat$ytil.obs <- pan.dat$y - pan.dat$x * coef(pols.obs)["x"]
  pols2.obs <- lm(ytil.obs ~ l2y + lx + lzo, data = pan.dat)


  ps.mod <- glm(x ~ ly + z + lx, data = pan.dat, na.action = na.exclude,
                family = binomial())
  num.mod <- glm(x ~  lx, data = pan.dat, na.action = na.exclude,
                 family = binomial())
  pscores <- fitted(ps.mod) * pan.dat$x + (1 - fitted(ps.mod)) * (1 - pan.dat$x)
  nscores <- fitted(num.mod) * pan.dat$x + (1 - fitted(num.mod)) * (1 - pan.dat$x)
  ws <- nscores/pscores
  cws <- pan.prod(ws, pan.dat$ii)
  msm <- lm(y ~ x + lx, data = pan.dat, weights = cws)

  raw.mod <- lm(y ~ x + lx, data = pan.dat)

  ps.obs.mod <- glm(x ~ ly + z.obs + lx, data = pan.dat, na.action = na.exclude,
                family = binomial())
  pscores.obs <- fitted(ps.obs.mod) * pan.dat$x + (1 - fitted(ps.obs.mod)) * (1 - pan.dat$x)
  ws.obs <- nscores/pscores.obs
  cws.obs <- pan.prod(ws.obs, pan.dat$ii)
  msm.obs <- lm(y ~ x + lx, data = pan.dat, weights = cws.obs)



  out <- rep(NA, times = 14)
  out[1] <- coef(pols)["lx"] + coef(pols)["ly"] * coef(pols)["x"]
  out[2] <- coef(pols2)["lx"]
  out[3] <- coef(msm)["lx"]
  out[4] <- coef(raw.mod)["lx"]
  out[5] <- mean(colMeans(y10-y00)[-c(1:2)])
  out[6] <- coef(pols)["x"]
  out[7] <- coef(msm)["x"]
  out[8] <- coef(raw.mod)["x"]
  out[9] <- mean(colMeans(y11-y10)[-c(1:2)])
  out[10] <- coef(pols.obs)["lx"] + coef(pols.obs)["ly"] * coef(pols.obs)["x"]
  out[11] <- coef(pols2.obs)["lx"]
  out[12] <- coef(msm.obs)["lx"]
  out[13] <- coef(pols.obs)["x"]
  out[14] <- coef(msm.obs)["x"]
  out[15] <- coef(pols.alt)["lx"] + coef(pols.alt)["ly"] * coef(pols.alt)["x"]

  if ((s %% 100) == 0) cat(s, "...")
  if ((s %% 1000) == 0) cat("\n")

  return(out)
}



## set up N/T parameters
Ns <- c(20, 30, 50, 100)
Ts <- c(10, 20, 50, 500)
params <- expand.grid(Ns, Ts)

## First
gamma1 <- c(0.5, 0)
gamma2 <- c(0.5, -0.5)
alpha1 <- c(-0.23, 2.5)
alpha2 <- c(-1.3, 2.5, 1.5)
mu1.1 <- 0
mu2.11 <- -0.1
mu2.01 <- -0.1
beta2.1 <- c(-0.1, 0)
beta2.0 <- c(-0.1, 0)
sig <- 0.1

sims <- 1000

## test
test.out <- tscs.sim(6, N = 20, T = 5, gamma = gamma1)
test.out

truth <- c(rep(mu1.1, 5), rep(mu2.11, 4), rep(mu1.1, 3), rep(mu2.11, 2), mu1.1)
meth.names <- c("ADL (t-1)", "SNMM (t-1)", "MSM (t-1)", "Raw (t-1)", "Truth (t-1)", "ADL/SNMM (t)", "MSM (t)", "Raw (t)", "Truth (t)", "ADL (t-1, mis)", "SNMM (t-1, mis)", "MSM (t-1, mis)", "ADL/SNMM (t, mis)", "MSM (t, mis)", "ADL (t-1, no Z)")
dims <- c(15, length(Ns), length(Ts))

bias.tvc.array <- array(dim = dims)
sds.tvc.array <- array(dim = dims)
rmse.tvc.array <- array(dim = dims)
bias.notvc.array <- array(dim = dims)
sds.notvc.array <- array(dim = dims)
rmse.notvc.array <- array(dim = dims)

dimnames(bias.tvc.array) <- list(meth.names, Ns, Ts)
dimnames(sds.tvc.array) <- list(meth.names, Ns, Ts)
dimnames(rmse.tvc.array) <- list(meth.names, Ns, Ts)
dimnames(bias.notvc.array) <- list(meth.names, Ns, Ts)
dimnames(sds.notvc.array) <- list(meth.names, Ns, Ts)
dimnames(rmse.notvc.array) <- list(meth.names, Ns, Ts)


for (i in 1:nrow(params)) {
  N.ind <- which(Ns == params[i,1])
  T.ind <- which(Ts == params[i,2])

  sim.out <- parallel::mclapply(1:sims, FUN = tscs.sim, mc.cores = 12,
                                N = params[i,1], T = params[i,2],
                                gamma = gamma1)
  sim.out.mat <- t(simplify2array(sim.out))
  bias <- colMeans(sim.out.mat) - truth
  sds <- apply(sim.out.mat, 2, sd)
  names(bias) <- meth.names
  names(sds) <- meth.names
  rmse <- sqrt(bias^2 + sds^2)
  bias.notvc.array[, N.ind, T.ind] <- bias
  sds.notvc.array[, N.ind, T.ind] <- sds
  rmse.notvc.array[, N.ind, T.ind] <- rmse
  cat("## No TVC, N=", params[i,1], " T=", params[i,2], " ##\n")

  print(cbind(bias=100*bias, sd=100*sds, rmse=100*rmse))

  sim.out <- parallel::mclapply(1:sims, FUN = tscs.sim, mc.cores = 12,
                                N = params[i,1], T = params[i,2],
                                gamma = gamma2)
  sim.out.mat <- t(simplify2array(sim.out))
  bias <- colMeans(sim.out.mat) - truth
  sds <- apply(sim.out.mat, 2, sd)
  names(bias) <- meth.names
  names(sds) <- meth.names
  rmse <- sqrt(bias^2 + sds^2)
  bias.tvc.array[, N.ind, T.ind] <- bias
  sds.tvc.array[, N.ind, T.ind] <- sds
  rmse.tvc.array[, N.ind, T.ind] <- rmse
  cat("## w/ TVC, N=", params[i,1], " T=", params[i,2], " ##\n")

  print(cbind(bias=100*bias, sd=100*sds, rmse=100*rmse))
}

cols <- brewer.pal(7, "Dark2")
cols <- c(cols[1:5], cols[2:5], cols[1:3], cols[2:3])
pchs <- c(15:18, 4, 16:18, 4, 15:17, 16:17)
ltys <- c(1:5, 2:5, 1:3, 2:3)


cairo_pdf(filename=here::here("output/fig5-sim.pdf"), family = "Work Sans Regular", width = 10, height = 10, pointsize = 14)
par(mfrow = c(2,2), oma = c(2, 0.5, 0.5, 0.5))
for (j in 2:3) {
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(rmse.tvc.array[1:5,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "RMSE", main = paste("Z endogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 1:5) {
    lines(Ns, abs(rmse.tvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(rmse.tvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(rmse.tvc.array[1:5,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "RMSE", main = paste("Z exogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 1:5) {
    lines(Ns, abs(rmse.notvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(rmse.notvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", meth.names[1:5], xpd = TRUE, horiz = TRUE, bty = "n", col = cols, pch = pchs, lty = 1:5)
dev.off()

cairo_pdf(filename=here::here("output/adl-snmm-msm-bias.pdf"), family = "Work Sans Regular", width = 10, height = 10, pointsize = 14)
par(mfrow = c(2,2), oma = c(2, 0.5, 0.5, 0.5))
for (j in 2:3) {
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(bias.tvc.array[1:5,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "Abs. Bias", main = paste("Z endogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 1:5) {
    lines(Ns, abs(bias.tvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(bias.tvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(bias.tvc.array[1:5,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "Abs. Bias", main = paste("Z exogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 1:5) {
    lines(Ns, abs(bias.notvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(bias.notvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", meth.names[1:5], xpd = TRUE, horiz = TRUE, bty = "n", col = cols[1:5], pch = pchs[1:5], lty = ltys[1:5])
dev.off()

cairo_pdf(filename=here::here("output/adl-snmm-msm-t-rmse.pdf"), family = "Work Sans Regular", width = 10, height = 5, pointsize = 14)
par(mfrow = c(1,2), oma = c(2, 0.5, 0.5, 0.5))
plot(NA, type = "o", xlim = c(10,500), ylim = range(abs(rmse.tvc.array[1:5,,])), las = 1,
       xlab = "T", ylab = "RMSE", main = paste("Z endogeneous (N = ", Ns[3], ")", sep = ""))
for (i in 1:5) {
  lines(Ts, abs(rmse.tvc.array[i,3,]), lty = ltys[i], col = cols[i])
  points(Ts, abs(rmse.tvc.array[i,3,]), pch = pchs[i], col = cols[i])
}
plot(NA, type = "o", xlim = c(10,500), ylim = range(abs(rmse.tvc.array[1:5,,])), las = 1,
     xlab = "T", ylab = "RMSE", main = paste("Z exogeneous (N = ", Ns[3], ")", sep = ""))
for (i in 1:5) {
  lines(Ts, abs(rmse.notvc.array[i,3,]), lty = ltys[i], col = cols[i])
  points(Ts, abs(rmse.notvc.array[i,3,]), pch = pchs[i], col = cols[i])
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", meth.names[1:5], xpd = TRUE, horiz = TRUE, bty = "n", col = cols, pch = pchs, lty = 1:5)
dev.off()

cairo_pdf(filename=here::here("output/adl-snmm-msm-t-bias.pdf"), family = "Work Sans Regular", width = 10, height = 5, pointsize = 14)
par(mfrow = c(1,2), oma = c(2, 0.5, 0.5, 0.5))
plot(NA, type = "o", xlim = c(10,500), ylim = range(abs(bias.tvc.array[1:5,,])), las = 1,
       xlab = "T", ylab = "Abs Bias", main = paste("Z endogeneous (N = ", Ns[3], ")", sep = ""))
for (i in 1:5) {
  lines(Ts, abs(bias.tvc.array[i,3,]), lty = ltys[i], col = cols[i])
  points(Ts, abs(bias.tvc.array[i,3,]), pch = pchs[i], col = cols[i])
}
plot(NA, type = "o", xlim = c(10,500), ylim = range(abs(bias.tvc.array[1:5,,])), las = 1,
     xlab = "T", ylab = "Abs Bias", main = paste("Z exogeneous (N = ", Ns[3], ")", sep = ""))
for (i in 1:5) {
  lines(Ts, abs(bias.notvc.array[i,3,]), lty = ltys[i], col = cols[i])
  points(Ts, abs(bias.notvc.array[i,3,]), pch = pchs[i], col = cols[i])
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", meth.names[1:5], xpd = TRUE, horiz = TRUE, bty = "n", col = cols, pch = pchs, lty = 1:5)
dev.off()



cairo_pdf(filename=here::here("output/adl-snmm-msm-cet-rmse.pdf"), family = "Work Sans Regular", width = 10, height = 10, pointsize = 14)
par(mfrow = c(2,2), oma = c(2, 0.5, 0.5, 0.5))
for (j in 2:3) {
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(rmse.tvc.array[6:9,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "RMSE", main = paste("Z endogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 6:9) {
    lines(Ns, abs(rmse.tvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(rmse.tvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(rmse.tvc.array[6:9,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "RMSE", main = paste("Z exogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 6:9) {
    lines(Ns, abs(rmse.notvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(rmse.notvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", meth.names[6:9], xpd = TRUE, horiz = TRUE, bty = "n", col = cols[6:9], pch = pchs[6:9], lty = ltys[6:9])
dev.off()

cairo_pdf(filename=here::here("output/adl-snmm-msm-cet-bias.pdf"), family = "Work Sans Regular", width = 10, height = 10, pointsize = 14)
par(mfrow = c(2,2), oma = c(2, 0.5, 0.5, 0.5))
for (j in 2:3) {
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(bias.tvc.array[6:9,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "Abs. Bias", main = paste("Z endogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 6:9) {
    lines(Ns, abs(bias.tvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(bias.tvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(bias.tvc.array[6:9,,])), las = 1,
       xlab = "Sample Size", ylab = "Abs. Bias", main = paste("Z exogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in 6:9) {
    lines(Ns, abs(bias.notvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(bias.notvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", meth.names[6:9], xpd = TRUE, horiz = TRUE, bty = "n", col = cols[6:9], pch = pchs[6:9], lty = ltys[6:9])
dev.off()



mis.lag.ind <- c(10,11,12,4,5)
par(mfrow = c(2,2), oma = c(2, 0.5, 0.5, 0.5))
for (j in 2:3) {
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(rmse.tvc.array[mis.lag.ind,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "RMSE", main = paste("Z endogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in mis.lag.ind) {
    lines(Ns, abs(rmse.tvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(rmse.tvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(rmse.tvc.array[mis.lag.ind,,])), las = 1,
       xlab = "Sample Size (N)", ylab = "RMSE", main = paste("Z exogeneous (T = ", Ts[j], ")", sep = ""))
  for (i in mis.lag.ind) {
    lines(Ns, abs(rmse.notvc.array[i,,j]), lty = ltys[i], col = cols[i])
    points(Ns, abs(rmse.notvc.array[i,,j]), pch = pchs[i], col = cols[i])
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", meth.names[mis.lag.ind], xpd = TRUE, horiz = TRUE, bty = "n", col = cols[mis.lag.ind], pch = pchs[mis.lag.ind], lty = ltys[mis.lag.ind])

ind.list <- list(lag.corr = 1:5,
                 lag.miss = c(10,11,12,4,5),
                 cet.corr = 6:9,
                 cet.miss = c(13,14, 8, 9))
ind.labs <- c("Lagged effect, correct Z",
              "Lagged effect, misspecified Z",
              "Contemporaneous effect, correct Z",
              "Contemporaneous effect, misspecified Z")
cairo_pdf(filename=here::here("output/miss-z-rmse.pdf"), family = "Work Sans Regular", width = 10, height = 10, pointsize = 14)
par(mfrow = c(2,2), oma = c(2, 0.5, 0.5, 0.5))
for (j in 1:4) {
  inds <- ind.list[[j]]
  plot(NA, type = "o", xlim = c(20,100), ylim = range(abs(rmse.tvc.array)), las = 1,
       xlab = "Sample Size (N)", ylab = "RMSE", main = ind.labs[[j]])
  for (i in inds) {
    lines(Ns, abs(rmse.tvc.array[i,,3]), lty = ltys[i], col = cols[i])
    points(Ns, abs(rmse.tvc.array[i,,3]), pch = pchs[i], col = cols[i])
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("ADL", "SNMM", "MSM", "Raw", "Truth"), xpd = TRUE, horiz = TRUE, bty = "n", col = cols[1:5], pch = pchs[1:5], lty = ltys[1:5])
dev.off()
