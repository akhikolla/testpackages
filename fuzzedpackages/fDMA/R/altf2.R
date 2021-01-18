
altf2 <- function (y,x,mods.incl=NULL,gprob=NULL,omega=NULL,av=NULL,window=NULL,initial.period=NULL,d=NULL,f=NULL,fmod=NULL,parallel=NULL)
  {

### computes some forecast quality measures for some alternative forecasts
### similar to altf()

### requires "forecast", "parallel", "stats" and "xts" packages

### y - a numeric or a column matrix of a dependent variable,

### x - a matrix of independent variables (drivers), different columns correspond to different independent variables

### mods.incl - a matrix indicating which models are used in averaging,
###             similar as in fDMA(),
###             if not specified, all possible models are used

### gprob - a matrix of Google probabilities, columns should correspond to columns of x

### omega - a parameter between 0 and 1 used in probabilities estimations,
###         used if gprob is specified

### av - models averaging method:
###      "ord" - each model is given the same weight,
###      "aic" - information-theoretic model averaging based on Akaike Information Criterion,
###      "aicc" - information-theoretic model averaging based on Akaike information Criterion with a correction for finite sample sizes,
###      "bic" - model averaging based on Bayesian Information Criterion,
###      "mse" - weights are computed according to Mean Squared Error (MSE)
###      if not specified, by default "ord" is used

### window - a size of a rolling regression window, a number of observations,
###          if not specified 10% of all observations are taken

### initial.period - a number of observation since which forecast quality measures are computed,
###                  by default the whole sample is used, i.e., initial.period = 1

### d - logical, used for hit.ratio calculation,
###     d = FALSE for level time-series,
###     d = TRUE if time-series represent changes,
###     by default d = FALSE

### f - vector of logical arguments indicating which forecast will be computed

### fmod - estimated model, class "dma" object

### parallel - indicate whether parallel computations should be used,
###            by default parallel = FALSE


### checking initial data

if (missing(y)) { stop("please, specify y") }
if (missing(x)) { stop("please, specify x") }
if (! (is.numeric(y) || is.matrix(y))) { stop("y must be numeric or matrix") }
if (is.matrix(y) && ! (ncol(y) == 1)) { stop("y must be a one column matrix") }
if (! is.matrix(x)) { stop("x must be a matrix") }
if (is.null(colnames(x)))
  {
    colnames(x) <- colnames(x, do.NULL = FALSE, prefix = "X")
    warning('column names of x were automatically created')
  }
if (anyNA(colnames(x))) { stop("x must have column names") }
if (is.matrix(y) && is.null(colnames(y)))
  {
    warning('column name of y was automatically created')
    colnames(y) <- colnames(y, do.NULL = FALSE, prefix = "Y")
  }
if (is.matrix(y) && anyNA(colnames(y)))
  {
    warning('column name of y was automatically created')
    colnames(y) <- "Y1"
  }
if (! length(y) == nrow(x)) { stop("y and x must have the same number of observations") }
if (anyNA(y)) { stop("missing values in y") }
if (anyNA(x)) { stop("missing values in x") }
if (is.null(mods.incl))
  {
    mods.incl <- expand.grid(rep.int(list(0:1), ncol(x)))
    mods.incl <- as.matrix(cbind(rep.int(1,nrow(mods.incl)),mods.incl))
  }
if ((! is.null(mods.incl)) && (! is.matrix(mods.incl))) { stop("mods.incl must be a matrix") }
if (is.matrix(mods.incl) && (! (ncol(mods.incl) == (ncol(x)+1)))) { stop("columns of mods.incl do not correspond to variables specified by columns of x") }
if ( is.matrix(mods.incl) && (length(mods.incl[!(mods.incl %in% c(0,1))]) > 0) ) { stop("mods.incl should contain only 0 and 1") }
if (is.matrix(mods.incl) && any(duplicated(mods.incl))) { stop("mods.incl contain duplicated models") }
if (is.matrix(mods.incl))
  {
    test <- FALSE
    test.row <- rep.int(0,ncol(mods.incl))
    for (i in 1:nrow(mods.incl))
      {
        if (identical(test.row,mods.incl[i,])) { test <- TRUE }
      }
    if (test == TRUE) { stop("mods.incl contain a model with no variables") }
  }
if (is.null(window)) { window <- floor(length(y)/10) }
if (window < 1) { window <- 1 }
if (!(is.numeric(window))) { stop("window must be numeric") }
if ((window < 0) || (window > length(y))) { stop("window must be a positive number less then the total number of observations") }
if (length(y) < 3) { stop("time-series too short: there have to be more than 3 observations") }
if (is.null(initial.period)) { initial.period <- 1 }
if (! is.numeric(initial.period)) { stop("initial.period must be numeric") }
if ((initial.period <= 0) || (initial.period > length(y))) { stop("initial.period must be greater than or equal to 1, and less than the number of observations") }
if (is.null(d)) { d <- FALSE }
if (! is.logical(d)) { stop("d must be logical, i.e., TRUE or FALSE") }
if (is.null(f)) { f <- rep(TRUE,4) }
if (is.null(av)) { av <- "ord" }
if (! av %in% c("ord","aic","aicc","bic","mse")) { stop("please, specify correct models averaging method") }
if (is.null(parallel)) { parallel <- FALSE }
if (! is.logical(parallel)) { stop("parallel must be logical, i.e., TRUE or FALSE") }

y <- as.matrix(y)
x <- as.matrix(x)

if (parallel == TRUE)
  {
     cl <- makeCluster(detectCores() - 1)
     clusterEvalQ(cl, {library(xts)})
     clusterExport(cl, c("y","x","mods.incl","av","window"), envir=environment())
  }

f.c <- function(ics)
  {
    for (i in 1:nrow(ics))
      {
        ics[i,] <- ics[i,] - min(ics[i,])
        ics[i,] <- exp(-ics[i,] / 2)
        ics[i,] <- ics[i,] / sum(ics[i,])
      }
    return(ics)
  }

f.gprob <- function(i)
  {
    if (length(which(mods.incl[i,-1,drop=FALSE]==1))>0)
      {
        p1 <- gprob[,which(mods.incl[i,-1,drop=FALSE]==1),drop=FALSE]
        gp1 <- apply(p1, 1, prod)
        gp1[gp1==0] <- 0.001 * (1/(2^ncol(mods.incl)))
      }
    else
      {
        gp1 <- matrix(1,ncol=1,nrow=nrow(gprob))
      }
    if (length(which(mods.incl[i,-1,drop=FALSE]==0))>0)
      {
        p2 <- 1-gprob[,which(mods.incl[i,-1,drop=FALSE]==0),drop=FALSE]
        gp2 <- apply(p2, 1, prod)
        gp2[gp2==0] <- 0.001 * (1/(2^ncol(mods.incl)))
      }
    else
      {
        gp2 <- matrix(1,ncol=1,nrow=nrow(gprob))
      }

    return( gp1 * gp2 )
  }

gprob.old <- gprob
if (!is.null(gprob)) { gprob.old <- as.matrix(gprob.old) }

if (! is.null(gprob))
  {
    gprob.ind <- nrow(x) - nrow(gprob) + 1
  }
else
  {
    gprob.ind <- nrow(x) + 1
  }
if ( ! is.null(gprob) && !(length(y) == nrow(gprob)) ) { warning("time-series of gprob and x differ in length") }

######################### OLS
##################################################

if (f[1]==TRUE)
{

f.ols <- function(i)
  {
    ind <- which(mods.incl[i,-1,drop=FALSE]==1)
    c <- mods.incl[i,1]
    if (length(ind)>0)
      {
        x.i <- x[,ind,drop=FALSE]
        if (c==1)
          {
            m <- lm(y[1:(initial.period-1)] ~ x.i[1:(initial.period-1),,drop=FALSE])
            xx <- cbind(1,x.i)
          }
        else
          {
            m <- lm(y[1:(initial.period-1)] ~ x.i[1:(initial.period-1),,drop=FALSE] -1)
            xx <- x.i
          }
      }
    else
      {
        m <- lm(y[1:(initial.period-1)] ~ 1)
        xx <- matrix(1,ncol=1,nrow=length(y))
      }
    
    fv <- as.vector(t(as.matrix(m$coefficients)) %*% t(xx))    
    aic <- AIC(m)
    n.par <- length(ind) + c + 1
    aicc <- aic + (2*n.par*(n.par+1))/(length(y)-n.par-1)
    bic <- BIC(m)
    mse <- mean((m$residuals)^2)
    mm <- summary(m)
    if (all(is.finite(m$coefficients)) && all(is.finite(mm$coefficients[,4])))
      {
        if (c==1)
          {
            coeff.m <- mm$coefficients[,1]
            pval.m <- mm$coefficients[,4]
          }
        else
          {
            coeff.m <- c(0,mm$coefficients[,1])
            pval.m <- c(1,mm$coefficients[,4])
          }
      }
    else
      {
        coeff.m <- rep(0,length(ind)+1)
        pval.m <- rep(1,length(ind)+1)
      }
    coeff <- matrix(0,nrow=1,ncol=ncol(x)+1)
    pval <- matrix(1,nrow=1,ncol=ncol(x)+1)
    coeff[,1] <- coeff.m[1]
    pval[,1] <- pval.m[1]
    if (length(coeff.m)>1)
      {
        coeff[,1+ind] <- coeff.m[-1]
        pval[,1+ind] <- pval.m[-1]
      }
    return(list(fv,aic,aicc,bic,mse,coeff,pval))
  }

if (parallel == TRUE)
  {
    y.ols <- parLapply(cl,seq(nrow(mods.incl)),f.ols)
  }
else
  {
    y.ols <- lapply(seq(nrow(mods.incl)),f.ols)
  }

if (av == "ord")
  {
    w <- rep.int(1 / nrow(mods.incl),nrow(mods.incl))
  }
if (av == "aic")
  {
    w <- sapply(y.ols,"[[",2)
    w <- w - min(w)
    w <- exp(-w / 2)
    w <- w / sum(w)
  }
if (av == "aicc")
  {
    w <- sapply(y.ols,"[[",3)
    w <- w - min(w)
    w <- exp(-w / 2)
    w <- w / sum(w)
  }
if (av == "bic")
  {
    w <- sapply(y.ols,"[[",4)
    w <- w - min(w)
    w <- exp(-w / 2)
    w <- w / sum(w)
  }
if (av == "mse")
  {
    w <- sapply(y.ols,"[[",5)
    w <- 1 / w
    w <- w / sum(w)
  }

w <- as.matrix(w)
w[!is.finite(w)] <- NA

if (!is.null(gprob.old))
  {
    gprob <- gprob.old[1:(initial.period-1),,drop=FALSE]
    gprob <- t(as.matrix(colMeans(gprob)))
    w.g <- t(sapply(seq(nrow(mods.incl)),f.gprob))
    w.g <- t(gNormalize(w.g))
    w <- omega * w + (1-omega) * w.g
  }

coeff.ols <- t(sapply(y.ols,"[[",6))
pval.ols <- t(sapply(y.ols,"[[",7))
coeff.ols <- t(w) %*% coeff.ols
pval.ols <- t(w) %*% pval.ols
colnames(coeff.ols) <- c("const",colnames(x))
colnames(pval.ols) <- colnames(coeff.ols)

y.ols <- sapply(y.ols,"[[",1)
y.ols <- y.ols %*% w
y.ols <- as.vector(y.ols)

w.ols <- t(w)

post.inc.ols <- t(w) %*% mods.incl
colnames(post.inc.ols) <- colnames(coeff.ols)

exp.var.ols <- t(w) %*% (as.matrix(apply(mods.incl,1,sum)))

}
else
{
y.ols <- NULL
coeff.ols <- NULL
pval.ols <- NULL
w.ols <- NULL
post.inc.ols <- NULL
exp.var.ols <- NULL
}

######################### recursive OLS
##################################################

if (f[2]==TRUE)
{

f.rec.ols <- function(i)
  {
    coeff.all <- matrix(0,nrow=nrow(x),ncol=ncol(x)+1)
    pval.all <- matrix(1,nrow=nrow(x),ncol=ncol(x)+1)

    ind <- which(mods.incl[i,-1,drop=FALSE]==1)
    c <- mods.incl[i,1]
    x.i <- x[,ind,drop=FALSE]
    if (length(ind)>0)
      {
        if (c==1)
          {
            m <- rec.reg(y=y,x=x.i,c=TRUE)
          }
        else
          {
            m <- rec.reg(y=y,x=x.i,c=FALSE)
          }
      }
    else
      {
        m <- rec.reg(y=y)
      }

    if (c==1)
      {
        coeff.m <- m$coeff.
        pval.m <- m$p.val
      }
    else
      {
        coeff.m <- cbind(0,m$coeff.)
        pval.m <- cbind(1,m$p.val)
      }

    coeff.all[,1] <- coeff.m[,1]
    pval.all[,1] <- pval.m[,1]

    if (ncol(coeff.m)>1)
      {
        coeff.all[,1+ind] <- coeff.m[,-1]
        pval.all[,1+ind] <- pval.m[,-1]
      }

    return(list(m$y.hat,m$AIC,m$AICc,m$BIC,m$MSE,coeff.all,pval.all))
  }

if (parallel == TRUE)
  {
    y.rec.ols <- parLapply(cl,seq(nrow(mods.incl)),f.rec.ols)
  }
else
  {
    y.rec.ols <- lapply(seq(nrow(mods.incl)),f.rec.ols)
  }

if (av == "ord")
  {
    w <- matrix(1 / nrow(mods.incl),nrow=nrow(x),ncol=nrow(mods.incl))
  }
if (av == "aic")
  {
    w <- sapply(y.rec.ols,"[[",2)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "aicc")
  {
    w <- sapply(y.rec.ols,"[[",3)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "bic")
  {
    w <- sapply(y.rec.ols,"[[",4)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "mse")
  {
    w <- sapply(y.rec.ols,"[[",5)
    w <- 1 / w
    w <- gNormalize(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }

w <- as.matrix(w)
w[!is.finite(w)] <- NA

if (!is.null(gprob.old))
  {
    gprob <- t(sapply(seq(1,nrow(gprob.old)), function(i) colMeans(gprob.old[1:i,,drop=FALSE])))
    w.g <- sapply(seq(nrow(mods.incl)),f.gprob)
    w.g <- gNormalize(w.g)
    w[gprob.ind:nrow(w),] <- omega * w[gprob.ind:nrow(w),] + (1-omega) * w.g
  }

w.rec.ols <- w

post.inc.rec.ols <- w %*% mods.incl

exp.var.rec.ols <- w %*% (as.matrix(apply(mods.incl,1,sum)))

coeff <- lapply(y.rec.ols,"[[",6)
pval <- lapply(y.rec.ols,"[[",7)

f.coeff <- function(i)
  {
    return(coeff[[i]][t,])
  }

f.pval <- function(i)
  {
    return(pval[[i]][t,])
  }

coeff.rec.ols <- matrix(0,nrow=1,ncol=ncol(x)+1)
pval.rec.ols <- matrix(1,nrow=1,ncol=ncol(x)+1)

for (t in 1:nrow(x))
  {
    coeff.av <- t(sapply(seq(length(coeff)),f.coeff))
    pval.av <- t(sapply(seq(length(pval)),f.pval))
    coeff.av <- w[t,] %*% coeff.av
    pval.av <- w[t,] %*% pval.av
    coeff.rec.ols <- rbind(coeff.rec.ols,coeff.av)
    pval.rec.ols <- rbind(pval.rec.ols,pval.av)
  }
coeff.rec.ols <- coeff.rec.ols[-1,]
pval.rec.ols <- pval.rec.ols[-1,]

colnames(coeff.rec.ols) <- c("const",colnames(x))
colnames(pval.rec.ols) <- colnames(coeff.rec.ols)
colnames(post.inc.rec.ols) <- colnames(coeff.rec.ols)

y.rec.ols <- sapply(y.rec.ols,"[[",1)
y.rec.ols <- y.rec.ols * w
y.rec.ols <- rowSums(y.rec.ols)
y.rec.ols <- as.vector(y.rec.ols)

}
else
{
y.rec.ols <- NULL
coeff.rec.ols <- NULL
pval.rec.ols <- NULL
w.rec.ols <- NULL
post.inc.rec.ols <- NULL
exp.var.rec.ols <- NULL
}

######################### rolling OLS
##################################################

if (f[3]==TRUE)
{

f.roll.ols <- function(i)
  {
    ind <- which(mods.incl[i,-1,drop=FALSE]==1)
    c <- mods.incl[i,1]
    if (length(ind)>0)
      {
        x.i <- x[,ind,drop=FALSE]
        if (c==1)
          {
            m <- roll.reg(y=y,x=x.i,window=window,c=TRUE)
          }
        else
          {
            m <- roll.reg(y=y,x=x.i,window=window,c=FALSE)
          }
      }
    else
      {
        m <- roll.reg(y=y,x=NULL,window=window)
      }

    coeff.all <- matrix(0,nrow=nrow(x),ncol=ncol(x)+1)
    pval.all <- matrix(1,nrow=nrow(x),ncol=ncol(x)+1)

    if (c==1)
      {
        coeff.m <- m$coeff.
        pval.m <- m$p.val
      }
    else
      {
        coeff.m <- cbind(0,m$coeff.)
        pval.m <- cbind(1,m$p.val)
      }

    coeff.all[,1] <- coeff.m[,1]
    pval.all[,1] <- pval.m[,1]
    if (ncol(coeff.m)>1)
      {
        coeff.all[,1+ind] <- coeff.m[,-1]
        pval.all[,1+ind] <- pval.m[,-1]
      }

    return(list(m$y.hat,m$AIC,m$AICc,m$BIC,m$MSE,coeff.all,pval.all))
  }

if (parallel == TRUE)
  {
    y.roll.ols <- parLapply(cl,seq(nrow(mods.incl)),f.roll.ols)
  }
else
  {
    y.roll.ols <- lapply(seq(nrow(mods.incl)),f.roll.ols)
  }

if (av == "ord")
  {
    w <- matrix(1 / nrow(mods.incl),nrow=nrow(x),ncol=nrow(mods.incl))
  }
if (av == "aic")
  {
    w <- sapply(y.roll.ols,"[[",2)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "aicc")
  {
    w <- sapply(y.roll.ols,"[[",3)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "bic")
  {
    w <- sapply(y.roll.ols,"[[",4)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "mse")
  {
    w <- sapply(y.roll.ols,"[[",5)
    w <- 1 / w
    w <- gNormalize(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }

w <- as.matrix(w)
w[!is.finite(w)] <- NA

if (!is.null(gprob.old))
  {
    gprob1 <- t(sapply(seq(1,(window-1)), function(i) colMeans(gprob.old[1:i,,drop=FALSE])))
    gprob2 <- t(sapply(seq(window,nrow(gprob.old)), function(i) colMeans(gprob.old[(i-window+1):i,,drop=FALSE])))
    gprob <- rbind(gprob1,gprob2)
    w.g <- sapply(seq(nrow(mods.incl)),f.gprob)
    w.g <- gNormalize(w.g)
    w[gprob.ind:nrow(w),] <- omega * w[gprob.ind:nrow(w),] + (1-omega) * w.g
  }

w.roll.ols <- w

post.inc.roll.ols <- w %*% mods.incl

exp.var.roll.ols <- w %*% (as.matrix(apply(mods.incl,1,sum)))

coeff <- lapply(y.roll.ols,"[[",6)
pval <- lapply(y.roll.ols,"[[",7)

f.coeff <- function(i)
  {
    return(coeff[[i]][t,])
  }

f.pval <- function(i)
  {
    return(pval[[i]][t,])
  }

coeff.roll.ols <- matrix(0,nrow=1,ncol=ncol(x)+1)
pval.roll.ols <- matrix(1,nrow=1,ncol=ncol(x)+1)

for (t in 1:nrow(x))
  {
    coeff.av <- t(sapply(seq(length(coeff)),f.coeff))
    pval.av <- t(sapply(seq(length(pval)),f.pval))
    coeff.av <- w[t,] %*% coeff.av
    pval.av <- w[t,] %*% pval.av
    coeff.roll.ols <- rbind(coeff.roll.ols,coeff.av)
    pval.roll.ols <- rbind(pval.roll.ols,pval.av)
  }
coeff.roll.ols <- coeff.roll.ols[-1,]
pval.roll.ols <- pval.roll.ols[-1,]

colnames(coeff.roll.ols) <- c("const",colnames(x))
colnames(pval.roll.ols) <- colnames(coeff.roll.ols)
colnames(post.inc.roll.ols) <- colnames(coeff.roll.ols)


y.roll.ols <- sapply(y.roll.ols,"[[",1)
y.roll.ols <- y.roll.ols * w
y.roll.ols <- rowSums(y.roll.ols)
y.roll.ols <- as.vector(y.roll.ols)

}
else
{
y.roll.ols <- NULL
coeff.roll.ols <- NULL
pval.roll.ols <- NULL
w.roll.ols <- NULL
post.inc.roll.ols <- NULL
exp.var.roll.ols <- NULL
}

######################### TVP
##################################################

if (f[4]==TRUE)
{

f.tvp <- function(i)
  {
    ind <- which(mods.incl[i,-1,drop=FALSE]==1)
    c <- mods.incl[i,1]

    x.i <- x[,ind,drop=FALSE]

    if (c==1)
      {
        m <- tvp(y=y,x=x.i,V=1,lambda=0.99,c=TRUE)
        n.par <- ncol(x.i) + 2
      }
    else
      {
        m <- tvp(y=y,x=x.i,V=1,lambda=0.99,c=FALSE)
        n.par <- ncol(x.i) + 1
      }
    aic <- vector()
    aicc <- vector()
    bic <- vector()
    for (i in 1:length(as.vector(y)))
      {
        if (c==1)
          {
            aic[i] <- -2 * log(m$pred.dens.[i]) + 2 * (ncol(x.i) + 1)
            aicc[i] <- aic[i] + (2*n.par*(n.par+1))/(i-n.par-1)
            bic[i] <- -2 * log(m$pred.dens.[i]) + (ncol(x.i) + 1) * log(i)
          }
        else
          {
            aic[i] <- -2 * log(m$pred.dens.[i]) + 2 * ncol(x.i)
            aicc[i] <- aic[i] + (2*n.par*(n.par+1))/(i-n.par-1)
            bic[i] <- -2 * log(m$pred.dens.[i]) + ncol(x.i) * log(i)
          }
      }
    mse <- vector()
    fv <- m$y.hat
    mse[1] <- (fv[1]-as.vector(y)[1])^2
    for (j in 2:length(fv))
      {
        mse[j] <- mean((fv[1:j]-as.vector(y)[1:j])^2)
      }

    coeff.all <- matrix(0,nrow=nrow(x),ncol=ncol(x)+1)

    if (c==1)
      {
        coeff.m <- m$thetas
      }
    else
      {
        coeff.m <- cbind(0,m$thetas)
      }

    coeff.all[,1] <- coeff.m[,1]
    if (ncol(coeff.m)>1)
      {
        coeff.all[,1+ind] <- coeff.m[,-1]
      }

    return(list(fv,aic,aicc,bic,mse,coeff.all))
  }

if (parallel == TRUE)
  {
    y.tvp <- parLapply(cl,seq(nrow(mods.incl)),f.tvp)
  }
else
  {
    y.tvp <- lapply(seq(nrow(mods.incl)),f.tvp)
  }

if (av == "ord")
  {
    w <- matrix(1 / nrow(mods.incl),nrow=nrow(x),ncol=nrow(mods.incl))
  }
if (av == "aic")
  {
    w <- sapply(y.tvp,"[[",2)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "aicc")
  {
    w <- sapply(y.tvp,"[[",3)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "bic")
  {
    w <- sapply(y.tvp,"[[",4)
    w <- f.c(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }
if (av == "mse")
  {
    w <- sapply(y.tvp,"[[",5)
    w <- 1 / w
    w <- gNormalize(w)
    w <- rbind(rep.int(1 / nrow(mods.incl),nrow(mods.incl)), w)
    w <- w[-nrow(w),]
  }

w <- as.matrix(w)
w[!is.finite(w)] <- NA

if (!is.null(gprob.old))
  {
    gprob <- gprob.old
    w.g <- sapply(seq(nrow(mods.incl)),f.gprob)
    w.g <- gNormalize(w.g)
    w[gprob.ind:nrow(w),] <- omega * w[gprob.ind:nrow(w),] + (1-omega) * w.g
  }

w.tvp <- w

post.inc.tvp <- w %*% mods.incl

exp.var.tvp <- w %*% (as.matrix(apply(mods.incl,1,sum)))

coeff <- lapply(y.tvp,"[[",6)

f.coeff <- function(i)
  {
    return(coeff[[i]][t,])
  }

coeff.tvp <- matrix(0,nrow=1,ncol=ncol(x)+1)

for (t in 1:nrow(x))
  {
    coeff.av <- t(sapply(seq(length(coeff)),f.coeff))
    coeff.av <- w[t,] %*% coeff.av
    coeff.tvp <- rbind(coeff.tvp,coeff.av)
  }
coeff.tvp <- coeff.tvp[-1,]

colnames(coeff.tvp) <- c("const",colnames(x))
colnames(post.inc.tvp) <- colnames(coeff.tvp)

pval.tvp <- NA

y.tvp <- sapply(y.tvp,"[[",1)
y.tvp <- y.tvp * w
y.tvp <- rowSums(y.tvp)
y.tvp <- as.vector(y.tvp)

}
else
{
y.tvp <- NULL
coeff.tvp <- NULL
pval.tvp <- NULL
w.tvp <- NULL
post.inc.tvp <- NULL
exp.var.tvp <- NULL
}

##################################################



fq <- list(y.ols,
           y.rec.ols,
           y.roll.ols,
           y.tvp
           )

coeff <- list(coeff.ols,
              coeff.rec.ols,
              coeff.roll.ols,
              coeff.tvp
              )

pval <- list(pval.ols,
             pval.rec.ols,
             pval.roll.ols,
             pval.tvp
             )

weights <- list(w.ols,
                w.rec.ols,
                w.roll.ols,
                w.tvp
                )

post.inc <- list(post.inc.ols,
                 post.inc.rec.ols,
                 post.inc.roll.ols,
                 post.inc.tvp
                 )

exp.var <- list(exp.var.ols,
                exp.var.rec.ols,
                exp.var.roll.ols,
                exp.var.tvp
                )

fq2 <- fq[f]

for (i in 1:4)
{
  if (!is.null(fq[[i]]))
    {
      fq[[i]] <- c(
                   as.numeric(accuracy(f=(as.vector(fq[[i]]))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))])),
                   as.numeric(hit.ratio(y=as.vector(y),y.hat=as.vector(fq[[i]]),d=d))
                  )
    }
}

fq <- fq[f]

fq <- matrix(unlist(fq),ncol=6,byrow=TRUE)

rnames <- c("av. OLS","av. rec. OLS","av. roll. OLS","av. TVP")
rownames(fq) <- rnames[f]
names(fq2) <- rnames[f]
colnames(fq) <- c("ME","RMSE","MAE","MPE","MAPE","HR")

coeff <- coeff[f]
pval <- pval[f]
weights <- weights[f]
post.inc <- post.inc[f]
exp.var <- exp.var[f]
names(coeff) <- rnames[f]
names(pval) <- rnames[f]
names(weights) <- rnames[f]
names(post.inc) <- rnames[f]
names(exp.var) <- rnames[f]

if (! is.null(fmod))
  {
    y.dma <- fmod$y.hat
    a.dma <- c(
               as.numeric(accuracy(f=(as.vector(y.dma))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))])),
               as.numeric(hit.ratio(y=as.vector(y),y.hat=as.vector(y.dma),d=d))
              )
    fq <- rbind(a.dma,fq)
    rownames(fq)[1] <- "est. model"
  }

if (parallel == TRUE)
  {
    stopCluster(cl)
    rm(cl)
  }

r <- list(round(fq,digits=4),fq2,as.matrix(y),coeff,weights,pval,post.inc,exp.var)
names(r) <- c("summary","y.hat","y","coeff.","weights","p.val.","rel.var.imp.","exp.var.")
class(r) <- "altf2"
return(r)

  }
