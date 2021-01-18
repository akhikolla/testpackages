
altf4 <- function (y,x,windows,V=NULL,alpha=NULL,lambda=NULL,initial.period=NULL,d=NULL,fmod=NULL,parallel=NULL,c=NULL,small.c=NULL)
  {

### computes some forecast quality measures for an alternative forecast,
### similar to altf() and fDMA(),
### i.e., the averaging like in fDMA() is performed over time-varying parameters rolling regressions with different windows sizes

### requires "forecast", "parallel", "stats" and "xts" packages

### y - a numeric or a column matrix of a dependent variable,

### x - a matrix of independent variables (drivers), different columns correspond to different independent variables

### windows - a vector of windows used in rolling regressions

### V - initial variance in the state space equation for the recursive moment estimator updating method,
###     as in the paper by Raftery et al. (2010),
###     if not specified V = 1 is used

### alpha - a forgetting factor between 0 and 1 used in probabilities estimations,
###         if not specified, lambda = 0.99 is used

### lambda - a forgetting factor between 0 and 1 used in variance approximations,
###          if not specified, lambda = 0.99 is used

### initial.period - a number of observation since which forecast quality measures are computed,
###                  by default the whole sample is used, i.e., initial.period = 1

### d - logical, used for hit.ratio calculation,
###     d = FALSE for level time-series,
###     d = TRUE if time-series represent changes,
###     by default d = FALSE

### fmod - estimated model, class "dma" object

### parallel - indicate whether parallel computations should be used,
###            by default parallel = FALSE

### c - a parameter indicating whether constant is included

### small.c - a small constant added to posterior model probabilities,
###           to prevent possible reduction them to 0 due to computational issues,
###           by default small.c is taken as in small constant as in Eq. (17)
###           in Raftery et al. (2010)

### checking initial data

if (missing(y)) { stop("please, specify y") }
if (! (is.numeric(y) || is.matrix(y))) { stop("y must be numeric or matrix") }
if (is.matrix(y) && ! (ncol(y) == 1)) { stop("y must be a one column matrix") }
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
if (anyNA(y)) { stop("missing values in y") }
if (length(y) < 3) { stop("time-series too short: there have to be more than 3 observations") }
if (is.null(V)) { V <- 1 }
if (is.null(lambda)) { lambda <- 0.99 }
if (is.null(alpha)) { alpha <- 0.99 }
if (is.null(initial.period)) { initial.period <- 1 }
if (! is.numeric(initial.period)) { stop("initial.period must be numeric") }
if ((initial.period <= 0) || (initial.period > length(y))) { stop("initial.period must be greater than or equal to 1, and less than the number of observations") }
if (is.null(d)) { d <- FALSE }
if (! is.logical(d)) { stop("d must be logical, i.e., TRUE or FALSE") }
if (is.null(parallel)) { parallel <- FALSE }
if (! is.logical(parallel)) { stop("parallel must be logical, i.e., TRUE or FALSE") }
if (! is.null(small.c) && ! is.numeric(small.c)) { stop("small.c must be a (small) number") }
if (! is.null(small.c) && (small.c<0)) { stop("small.c must be positive") }


y <- as.matrix(y)
x <- as.matrix(x)

if ( is.null(c) ) { c <- TRUE }
if ( ncol(x) == 0 ) { c <- TRUE }

windows <- unique(windows)
windows <- sort(windows, decreasing=FALSE)

if (parallel == TRUE)
  {
     cl <- makeCluster(detectCores() - 1)
     clusterEvalQ(cl, {library(xts)})
     clusterExport(cl, c("y","x","windows","V","lambda"), envir=environment())
  }


######################### rolling TVP
##################################################

f.roll.ols <- function(i)
  {
    window <- windows[i]
    pd <- vector()
    y.roll.ols <- vector()
    if (c==TRUE)
      {
        coeff <- matrix(NA,nrow=1,ncol=ncol(x)+1)
      }
    else
      {
        coeff <- matrix(NA,nrow=1,ncol=ncol(x))
      }
    for (i in 1:(window-1))
      {
        if (i==1)
          {
            m <- tvp(y=y[1],x=x[1,,drop=FALSE],V=V,lambda=lambda,c=c)
            y.roll.ols[1] <- m$y.hat[1]
            pd[1] <- m$pred.dens.[1]
            coeff <- rbind(coeff,m$thetas[1,])
          }
        else
          {
            m <- tvp(y=as.vector(y)[1:i],x=x[1:i,,drop=FALSE],V=V,lambda=lambda,c=c)
            y.roll.ols[i] <- m$y.hat[i]
            pd[i] <- m$pred.dens.[i]
            coeff <- rbind(coeff,m$thetas[i,])
          }
      }

    for (i in window:nrow(x))
      {
        m <- tvp(y=as.vector(y)[(i-window+1):i],x=x[(i-window+1):i,,drop=FALSE],V=V,lambda=lambda,c=c)
        y.roll.ols[i] <- m$y.hat[window]
        pd[i] <- m$pred.dens.[window]
        coeff <- rbind(coeff,m$thetas[window,])
      }

    return(list(y.roll.ols,pd,coeff[-1,,drop=FALSE]))
  }

if (parallel == TRUE)
  {
    y.roll.ols <- parLapply(cl,seq(length(windows)),f.roll.ols)
  }
else
  {
    y.roll.ols <- lapply(seq(length(windows)),f.roll.ols)
  }

##################################################

w <- sapply(y.roll.ols,"[[",2)
coeff <- lapply(y.roll.ols,"[[",3)
y.roll.ols <- sapply(y.roll.ols,"[[",1)
pi1 <- as.vector(rep.int(1/length(windows),length(windows)))
if (is.null(small.c))
  {
    c2 <- 0.001 * (1/(2^length(windows)))
  }
else
  {
    c2 <- small.c
  }
y.pred <- vector()

if (c==TRUE)
  {
    coeff.av.all <- matrix(NA,nrow=1,ncol=ncol(x)+1)
  }
else
  {
    coeff.av.all <- matrix(NA,nrow=1,ncol=ncol(x))
  }
weights <- matrix(NA,nrow=1,ncol=length(windows))

f.thetas <- function(i)
  {
    return(coeff[[i]][t,])
  }

for (t in 1:nrow(x))
  {
    pi2 <- (pi1^alpha + c2) / (sum((pi1)^alpha + c2))
    y.pred[t] <- crossprod(pi2,y.roll.ols[t,])
    coeff.av <- t(sapply(seq(length(coeff)),f.thetas))
    if (ncol(x)==0) { coeff.av <- t(coeff.av) }
    coeff.av <- pi2 %*% coeff.av
    coeff.av.all <- rbind(coeff.av.all,coeff.av)
    weights <- rbind(weights,pi2)
    pi1 <- (pi2 * w[t,]) / as.numeric(crossprod(pi2,w[t,]))
  }

y.roll.ols <- y.pred

exp.win <- weights[-1,] %*% windows

##################################################

coeff.av.all <- list(coeff.av.all[-1,,drop=FALSE])
if (ncol(x)==0) { colnames(coeff.av.all[[1]]) <- "const" }
names(coeff.av.all) <- "av. roll. TVP"
weights <- list(weights[-1,])
names(weights) <- "av. roll. TVP"
exp.win <- list(exp.win)
names(exp.win) <- "av. roll. TVP"
rownames(weights[[1]]) <- rownames(coeff.av.all[[1]])
colnames(weights[[1]]) <- windows

fq2 <- list(y.roll.ols)
names(fq2) <- "av. roll. TVP"

fq <- c(as.numeric(accuracy(f=(as.vector(y.roll.ols))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))])),
        as.numeric(hit.ratio(y=as.vector(y),y.hat=as.vector(y.roll.ols),d=d))
       )

fq <- matrix(unlist(fq),ncol=6,byrow=TRUE)

rownames(fq) <- "av. roll. TVP"
colnames(fq) <- c("ME","RMSE","MAE","MPE","MAPE","HR")

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

r <- list(round(fq,digits=4),fq2,as.matrix(y),coeff.av.all,weights,exp.win)
names(r) <- c("summary","y.hat","y","coeff.","weights","exp.win.")
class(r) <- "altf4"
return(r)

  }
