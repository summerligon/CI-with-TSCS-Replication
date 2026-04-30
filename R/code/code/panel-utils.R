library(plyr)
library(reshape)
night <- rgb(26,31,30, max=255)
beach <- rgb(227,223,186, max=255)
red60 <- rgb(1,.4,.4)
tangyblue <- rgb(108,189,181, max=255)
purp <- rgb(181.5,145.5, 141.5, max=255)


green1 <- rgb(178,179,159, max = 255)
lightbrown <- rgb(205,140,82, max = 255)

green2 <- rgb(200,201,181, max = 255)
green3 <- rgb(222,223,197, max = 255)
orange <- rgb(61,66,60, max = 255)
aqua <- rgb(240,236,201, max = 255)

lower.95.bound <- function(x) quantile(x, prob = 0.025)
upper.95.bound <- function(x) quantile(x, prob = 0.975)
regTable <- each(mean, sd, lower.95.bound, upper.95.bound)

dt.to.week <- function(x) {
  x <- as.numeric(format(x, "%W"))

}

                                        #
                                        # Function: tslag (lagging a vector)
                                        #
tslag <- function(x, d=1) {
  x <- as.vector(x)
  n <- length(x)
  c(rep(NA,d),x)[1:n]
}

pastmin <- function(x) {
  xt <- x
  xt[!is.na(xt)] <- cummin(na.omit(x))
  return(xt)
}
pastmax <- function(x) {
  xt <- x
  xt[!is.na(xt)] <- cummax(na.omit(x))
  return(xt)
}

pastsum <- function(x) {
  xt <- x
  xt[!is.na(xt)] <- cumsum(na.omit(x))
  return(xt)
}

pan.lag <- function(x,ind,length = 1) {
  lag.one <- function(w) {
    adds <- rep(NA, times = length)
    drops <- seq(length(w), length(w)-(length-1), by = -1)
    return(c(adds, w[-drops]))
  }
  unlist(tapply(x,ind,lag.one))
}

pan.diff <- function(x,ind,length = 1) {
  return(x - pan.lag(x, ind, length))
}

d <- pan.diff

pan.sum <- function(x, ind) {
  unlist(tapply(x,ind, function(x) {
    xt <- x
    xt[!is.na(xt)] <- cumsum(na.omit(x))
    return(xt)
  }))
}

pan.prod <- function(x,ind) {
  unlist(tapply(x,ind, function(x) {
    xt <- x
    xt[!is.na(xt)] <- cumprod(na.omit(x))
    return(xt)
  }))
}

pan.mean <- function(x,ind) {
  unlist(tapply(x, ind, function(x) rep(mean(x, na.rm=TRUE), length(x))))
}

pan.first <- function(x,ind) {
  unlist(tapply(x, ind, function(x) rep(ifelse(any(x),which(x)[1],NA),length(x)) ))
}


pan.min <- function(x,ind) {
  unlist(tapply(x, ind, function(x) rep(min(x, na.rm=TRUE), length(x))))
}


pan.cummin <- function(x, ind) {
  unlist(tapply(x,ind, function(x) {
    xt <- x
    xt[!is.na(xt)] <- cummin(na.omit(x))
    return(xt)
  }))
}



## ones is an argument to manually set certain observations as having
## weight 1. This is needed due to common support or positivity
## problems. Basically, in certain ranges of the covariates, there are
## either random or structural zeros. For negative advertising, a
## strutural zero occurs when there are no ads. Random zeros generally
## occur in the extremes of the polls and number of ads
## covariates. Extremely uncompetitive races never see negativity and
## extreme competitive race (in the number of ads) almost always go
## negative.

iptw <- function(denominator, numerator, data, id, time, family, subset,
                 ones = NULL, bal.covs = NULL, eval.balance = TRUE) {

  require(mgcv)
  require(MASS)

  if ("family" %in% class(family)) {
    shortfam <- family$family
  } else {
    shortfam <- family
  }

  supported.families <- c("binomial", "gaussian", "ologit")
  if (!(shortfam %in% supported.families)) {
    stop("That family is not supported, silly.")
  }

  if (missing(ones)) {
    ones <- rep(FALSE, nrow(data))
  }

  ## via the subset.data.frame function
  if (missing(subset)) {
    r <- TRUE
  } else {
    e <- substitute(subset)
    r <- eval(e, data, parent.frame())
    if (!is.logical(r)) {
      stop("'subset' must evaluate to logical")
    }
    r <- r & !is.na(r)
    if (sum(r) == 0) {
      stop("no observations in the subset")
    }
  }

  odata <- data
  data <- data[r,]
  ones <- ones[r]
  idmat <- data[,c(id, time)]
  d.fit <- rep(NA, nrow(data))
  n.fit <- rep(NA, nrow(data))
  names(d.fit) <- rownames(data)
  names(n.fit) <- rownames(data)
  trvar <- all.vars(denominator)[1]
  tr <- data[, trvar]

  if (shortfam %in% c("binomial", "gaussian")) {
    dem <- gam(denominator, data = data[!ones,], family = family)
    nom <- gam(numerator, data = data[!ones,], family = family)
  }
  if (shortfam == "ologit") {
    dem <- polr(denominator, data = data[!ones,], Hess = TRUE)
    nom <- polr(numerator, data = data[!ones,], Hess = TRUE)
  }

  if (length(nom) > 1) {
    d.nms <- rownames(model.matrix(dem))
    n.nms <- rownames(model.matrix(nom))
    if (shortfam == "binomial") {
      d.fit[d.nms] <- dem$fitted.values
      n.fit[n.nms] <- nom$fitted.values

      d.pr.tr <- ifelse(tr == 1, d.fit, 1-d.fit)
      n.pr.tr <- ifelse(tr == 1, n.fit, 1-n.fit)
    }
    if (shortfam == "gaussian") {
      dem.dens <- dnorm(model.frame(dem)[,1], dem$fitted.values,
                        sqrt(summary(dem)$dispersion)[1])
      nom.dens <- dnorm(model.frame(nom)[,1], nom$fitted.values,
                        sqrt(summary(nom)$dispersion)[1])
      d.fit[d.nms] <- dem.dens
      n.fit[n.nms] <- nom.dens

      d.pr.tr <- d.fit
      n.pr.tr <- n.fit
    }
    if (shortfam == "ologit") {
      d.indic <- model.matrix( ~ dem$model[,1]-1)
      n.indic <- model.matrix( ~ nom$model[,1]-1)
      d.fit[d.nms] <- rowSums(d.indic * dem$fitted.values)
      n.fit[n.nms] <- rowSums(n.indic * nom$fitted.values)
      d.pr.tr <- d.fit
      n.pr.tr <- n.fit
    }
  } else {
    d.nms <- rownames(model.matrix(dem))

    d.fit[d.nms] <- dem$fitted.values
    n.fit[TRUE] <- nom
    d.pr.tr <- ifelse(tr == 1, d.fit, 1-d.fit)
    n.pr.tr <- ifelse(tr == 1, n.fit, 1-n.fit)
  }

  sw <- n.pr.tr/d.pr.tr
  names(sw) <- rownames(data)
  sw[ones] <- 1
  sw <- pan.prod(sw, idmat[,1])
  wtab <- data.frame(idmat[,1], idmat[,2], sw = sw, d.pscore = d.fit,
                     n.pscore = n.fit)
  colnames(wtab)[1:2] <- c(id, time)
  out <- list(sw = sw, wtab = wtab, d.mod = dem, n.mod = nom, d.pscore = d.fit,
              n.pscore = n.fit, ivar = id, tvar = time, subset = r, ones = ones)
  class(out) <- "iptw"
  if (shortfam == "binomial" && eval.balance) {
    bal.out <- balance.iptw(out, data = odata, add.covs = bal.covs)
    out$balance <- bal.out
  }
  return(out)
}

summary.iptw <- function(x) {
  print(x$balance[,c(3,6)])
}

plot.iptw <- function(x, logw = TRUE, ...) {
  if (logw) {
    x$wtab$sw <- log(x$wtab$sw)
  }
  colnames(x$wtab)[1:2] <- c("ivar", "tvar")
  tsum <- cast(na.omit(melt(x$wtab, id.variables = c("ivar", "tvar"),
                                measure=c("sw"))), tvar ~ ., summary)
  subtle.boxplot(tsum, ...)
  if (logw) {
    abline(h = 0, col = "grey")
  } else {
    abline(h = 1, col = "grey")
  }
  invisible(tsum)
}

merge.iptw <- function(x, data) {
  data$sw <- data$d.pscore <- data$n.pscore <- NULL
  data <- merge(data, x$wtab, by = c(x$ivar, x$tvar), all.x = TRUE)
  return(data)
}

sens.adj <- function(y, a, a.fits, alpha, confound.func = function(x,y) {x}) {
  prob.astar <- ifelse(a == 0, a.fits, 1 - a.fits)
  adjust <- prob.astar
  alpha.out <- confound.func(alpha, a)
  ionic <- alpha.out * prob.astar
#  ionic <- t(matrix(alpha.out, ncol = 1) %*% prob.astar)

  adjust <- pastsum(ionic)

  y.alpha <- y - adjust
  return(y.alpha)
}



subtle.boxplot <- function(sum, xlab = "", ylab = "") {
  plot(x = NULL, y = NULL, pch = 19, axes = FALSE,
       ylim = c(min(sum$Min), max(sum$Max)), xlim =
       range(sum[,1]), xlab = xlab, ylab = ylab)
  axis(side = 1)
  axis(side = 2, las = 2)

  segments(x0 = sum[,1], x1 = sum[,1],
           y0 = sum$Min, y1 = sum$Max, lwd = 1,
           col = rgb(0.25, 0.25, 0.25, alpha = 0.25))
  rect(xleft = sum[,1]-0.4, xright = sum[,1]+0.4,
       ybottom = sum$X1st.Qu., ytop = sum$X3rd.Qu.,
       col = rgb(0.75, 0.75, 0.75), border = FALSE)
  points(x = sum[,1], y = sum$Min, pch = 19,
         col = rgb(0.75, 0.75, 0.75))
  points(x = sum[,1], y = sum$Max, pch = 19,
         col = rgb(0.75, 0.75, 0.75))
  segments(x0 = sum[,1]-0.4, x1 = sum[,1] + 0.4,
           y0 = sum$Mean, y1 = sum$Mean, lwd = 3)
  invisible()
}


balance.iptw <- function(obj, data, add.covs = NULL) {
  require(survey)


  d.form <- interpret.gam(obj$d.mod$formula)$fake.formula
  d.vars <- all.vars(d.form[[3]])
  tr.var <- all.vars(d.form)[1]
  n.form <- interpret.gam(obj$n.mod$formula)$fake.formula
  n.vars <- all.vars(n.form[[3]])
  n.form <- as.character(n.form)[3]
  covs <- d.vars[!(d.vars %in% n.vars)]
  covs <- unique(c(covs, add.covs))

  data <- data[obj$subset, c(obj$ivar, obj$tvar, tr.var, covs, n.vars)]
  #data$ones <- obj$ones
  data$sw <- obj$sw

  id.form <- as.formula(paste("~",obj$ivar))

  causal.des <- svydesign(ids = id.form, weights = ~sw,
                          data = data[!obj$ones & !is.na(obj$sw),])

  all.bal <- panel.balance(treat = tr.var, bal.covs = covs, n.form = n.form,
                           design = causal.des)
  return(all.bal)

}

panel.balance <- function(treat, bal.covs, n.form, design, subset = NULL) {

  if (!is.null(subset)) {
    subset <- eval(parse(text = subset))
  }
  unweighted <- matrix(NA, nrow = length(bal.covs), ncol = 2)
  weighted <- unweighted

  for (i in seq(along.with=bal.covs)) {
    new.form <- as.formula(paste(bal.covs[i], "~", treat, "+", n.form))
    mc <- design$call
    mc$weights <- as.formula("~1")
    unw.design <- eval(mc, parent.frame())

    mc[[1]] <- as.name("lm")
    mc$formula <- new.form
    mc$subset <- subset
    unw <- svyglm(new.form, subset = subset, design = unw.design)
    unweighted[i, ] <- summary(unw)$coef[treat, 1:2]

    wei <- svyglm(new.form, subset = subset, design = design)
    weighted[i, ] <- summary(wei)$coef[treat, 1:2]
  }
  out <- cbind(unweighted, unweighted[,1]/unweighted[,2],
               weighted, weighted[,1]/weighted[,2])

  rownames(out) <- bal.covs
  colnames(out) <- c("Unweighed Coef.", "Unweighted SE", "USCoef",
                     "Weighted Coef.", "Weighted SE", "WSCoef")
  return(out)
}


## confint for gee LM object
confint.geelm <- function(object, parm, level = 0.95, ...) {
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm))
        parm <- pnames
    else if (is.numeric(parm))
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- stats:::format.perc(a, 3)
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
        pct))
    ses <- summary(object)$coefficients[parm,2]
    ci[] <- cf[parm] + ses %o% fac
    ci
}

l <- pan.lag

snmm.var <- function(mods, blip.vars, data, unit.var) {
  WW <- lapply(mods, model.matrix)
  cowlists<- lapply(WW, function(x) unique(data[rownames(x), unit.var]))
  u.inds <- sapply(cowlists, function(x) unique(data[,unit.var]) %in% x)
  unicow <- unique(data[,unit.var])[rowSums(u.inds) == length(cowlists)]
  W.tilde <- WW[-length(WW)]
  for (w in 1:length(W.tilde)) W.tilde[[w]][,-blip.vars] <- 0
  kvec <- sapply(WW, ncol)
  K <- sum(kvec)
  iii <- cumsum(kvec)
  ii <- cumsum(c(0,kvec[-length(kvec)])) + 1
  Ghat <- matrix(0, nrow = K, ncol = K)
  meat <- list()
  for (i in 1:length(kvec)) {
    Ghat[ii[i]:iii[i], ii[i]:iii[i]] <- t(WW[[i]]) %*% WW[[i]]
    meat[[i]] <- apply(estfun(mods[[i]]), 2, function(x) tapply(x, data[rownames(WW[[i]]), unit.var], sum))
    meat[[i]] <- meat[[i]][which(unique(data[rownames(WW[[i]]), unit.var]) %in% unicow),]
    for (j in i:length(kvec)) {
      if (i < j) {
        nms <- intersect(rownames(WW[[i]]), rownames(WW[[j]]))
        Ghat[ii[j]:iii[j], ii[i]:iii[i]] <- t(WW[[j]][nms,]) %*% W.tilde[[i]][nms,]
      }
    }
  }
  meat.big <- do.call(cbind, meat)
  M <- length(unicow)
  NN <- nrow(WW[[1]])
  dfc <- (M/(M-1))*((NN-1)/(NN-K))
  vcv <- dfc*solve(crossprod(Ghat)) %*% t(Ghat) %*% crossprod(meat.big) %*% Ghat %*% solve(crossprod(Ghat))
  return(vcv)
}


snmm.var.single <- function(mods, blip.vars, data) {
  WW <- lapply(mods, model.matrix)
  W.tilde <- WW[-length(WW)]
  for (w in 1:length(W.tilde)) W.tilde[[w]][,-blip.vars] <- 0
  kvec <- sapply(WW, ncol)
  K <- sum(kvec)
  iii <- cumsum(kvec)
  ii <- cumsum(c(0,kvec[-length(kvec)])) + 1
  Ghat <- matrix(0, nrow = K, ncol = K)
  meat <- list()
  for (i in 1:length(kvec)) {
    Ghat[ii[i]:iii[i], ii[i]:iii[i]] <- t(WW[[i]]) %*% WW[[i]]
    meat[[i]] <- estfun(mods[[i]])
    for (j in i:length(kvec)) {
      if (i < j) {
        nms <- intersect(rownames(WW[[i]]), rownames(WW[[j]]))
        Ghat[ii[j]:iii[j], ii[i]:iii[i]] <- t(WW[[j]][nms,]) %*% W.tilde[[i]][nms,]
      }
    }
  }
  meat.big <- do.call(cbind, meat)
  NN <- nrow(WW[[1]])
  dfc <- ((NN-1)/(NN-K))
  vcv <- dfc*solve(crossprod(Ghat)) %*% t(Ghat) %*% crossprod(meat.big) %*% Ghat %*% solve(crossprod(Ghat))
  return(vcv)
}
