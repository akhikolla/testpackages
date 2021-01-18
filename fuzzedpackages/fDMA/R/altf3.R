
altf3 <- function (y,x=NULL,windows,av=NULL,initial.period=NULL,d=NULL,fmod=NULL,parallel=NULL,c=NULL)
  {

### computes some forecast quality measures for an alternative forecast,
### similar to altf(),
### i.e., the averaging is performed over rolling regressions with different windows sizes
 
### requires "forecast", "parallel", "stats" and "xts" packages

### y - a numeric or a column matrix of a dependent variable,

### x - a matrix of independent variables (drivers), different columns correspond to different independent variables

### av - models averaging method:
###      "ord" - each model is given the same weight,
###      "aic" - information-theoretic model averaging based on Akaike Information Criterion,
###      "aicc" - information-theoretic model averaging based on Akaike information Criterion with a correction for finite sample sizes,
###      "bic" - model averaging based on Bayesian Information Criterion, 
###      "mse" - weights are computed according to Mean Squared Error (MSE),
###      if av is a number then weights are computed with respect to window size, 
###      if not specified, by default "ord" is used

### windows - a vector of windows used in rolling regressions

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
if (is.null(initial.period)) { initial.period <- 1 }
if (! is.numeric(initial.period)) { stop("initial.period must be numeric") }
if ((initial.period <= 0) || (initial.period > length(y))) { stop("initial.period must be greater than or equal to 1, and less than the number of observations") }
if (is.null(av)) { av <- "ord" }
if (! av %in% c("ord","aic","aicc","bic","mse") && ! is.numeric(av)) { stop("please, specify correct models averaging method") }
if (is.null(d)) { d <- FALSE }
if (! is.logical(d)) { stop("d must be logical, i.e., TRUE or FALSE") }
if (is.null(parallel)) { parallel <- FALSE }
if (! is.logical(parallel)) { stop("parallel must be logical, i.e., TRUE or FALSE") }


y <- as.matrix(y)
if (! is.null(x)) { x <- as.matrix(x) }
if (is.null(c)) { c <- TRUE }
if (is.null(x)) { c <- TRUE }

windows <- unique(windows)
windows <- sort(windows, decreasing=FALSE)

if (parallel == TRUE)
  {
     cl <- makeCluster(detectCores() - 1)
     clusterEvalQ(cl, {library(xts)})
     clusterExport(cl, c("y","x","windows","c"), envir=environment())
  }


######################### rolling OLS
##################################################

f.roll.ols <- function(i)
  {
    return(roll.reg(y=y,x=x,window=windows[i],c=c)) 
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
  
if (av == "ord") 
  { 
    w <- matrix(1 / length(windows),nrow=nrow(x),ncol=length(windows)) 
  }
if (av == "aic") 
  { 
    w <- sapply(y.roll.ols,"[[",2) 
    w <- f.c(w)
    w <- rbind(rep.int(1 / length(windows),length(windows)), w)
    w <- w[-nrow(w),]
  }
if (av == "aicc") 
  { 
    w <- sapply(y.roll.ols,"[[",3) 
    w <- f.c(w)
    w <- rbind(rep.int(1 / length(windows),length(windows)), w)
    w <- w[-nrow(w),]
  }
if (av == "bic") 
  { 
    w <- sapply(y.roll.ols,"[[",4) 
    w <- f.c(w)
    w <- rbind(rep.int(1 / length(windows),length(windows)), w)
    w <- w[-nrow(w),]
  }
if (av == "mse") 
  { 
    w <- sapply(y.roll.ols,"[[",5) 
    w <- 1 / w
    w <- gNormalize(w)
    w <- rbind(rep.int(1 / length(windows),length(windows)), w)
    w <- w[-nrow(w),]
  }
if (is.numeric(av))
  {
    w <- windows^av
    w <- w / sum (w)
    w <- matrix(w,nrow=length(as.vector(y)),ncol=length(w),byrow=TRUE)
  }

w <- as.matrix(w)

exp.win <- w %*% windows

coeff <- lapply(y.roll.ols,"[[",6)
pval <- lapply(y.roll.ols,"[[",7)

f.thetas <- function(i)
  {
    return(coeff[[i]][t,])
  }

f.pval <- function(i)
  {
    return(pval[[i]][t,])
  }

if (c==TRUE)
  {
    if (!is.null(x)) 
      { 
        coeff.av.all <- matrix(NA,nrow=1,ncol=ncol(x)+1)
      }
    else
      {
        coeff.av.all <- matrix(NA,nrow=1,ncol=1)
      }
  }
else
  {
    coeff.av.all <- matrix(NA,nrow=1,ncol=ncol(x))
  }
pval.av.all <- coeff.av.all

for (t in 1:length(as.vector(y)))
  {
    if (is.null(x)) 
      {
        coeff.av <- sapply(seq(length(coeff)),f.thetas)
        pval.av <- sapply(seq(length(coeff)),f.pval)
      }
    else
      {
        coeff.av <- t(sapply(seq(length(coeff)),f.thetas))
        pval.av <- t(sapply(seq(length(coeff)),f.pval))
      }
    coeff.av <- w[t,] %*% coeff.av
    pval.av <- w[t,] %*% pval.av
    coeff.av.all <- rbind(coeff.av.all,coeff.av)
    pval.av.all <- rbind(pval.av.all,pval.av)
  }
coeff.av.all <- coeff.av.all[-1,,drop=FALSE] 
pval.av.all <- pval.av.all[-1,,drop=FALSE] 
if (is.null(x))
  {
   colnames(coeff.av.all) <- "const"
   colnames(pval.av.all) <- colnames(coeff.av.all)
  }
 
y.roll.ols <- sapply(y.roll.ols,"[[",1)
y.roll.ols <- y.roll.ols * w
y.roll.ols <- rowSums(y.roll.ols)
y.roll.ols <- as.vector(y.roll.ols)

##################################################

coeff.av.all <- list(coeff.av.all)
names(coeff.av.all) <- "av. roll. TVP"
pval.av.all <- list(pval.av.all)
names(pval.av.all) <- "av. roll. TVP"
exp.win <- list(exp.win)
names(exp.win) <- "av. roll. TVP"
w <- list(w)
names(w) <- "av. roll. TVP"
colnames(w[[1]]) <- windows
rownames(w[[1]]) <- rownames(coeff.av.all[[1]])

fq2 <- list(y.roll.ols)
names(fq2) <- "av. roll. OLS"

fq <- c(as.numeric(accuracy(f=(as.vector(y.roll.ols))[initial.period:length(as.vector(y))],x=(as.vector(y))[initial.period:length(as.vector(y))])),
        as.numeric(hit.ratio(y=as.vector(y),y.hat=as.vector(y.roll.ols),d=d))
       )

fq <- matrix(unlist(fq),ncol=6,byrow=TRUE)

rownames(fq) <- "av. roll. OLS"
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

r <- list(round(fq,digits=4),fq2,as.matrix(y),coeff.av.all,w,pval.av.all,exp.win)
names(r) <- c("summary","y.hat","y","coeff.","weights","p.val.","exp.win.")
class(r) <- "altf3"
return(r)

  }
  