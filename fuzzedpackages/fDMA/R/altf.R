
altf <- function (y,x,window=NULL,initial.period=NULL,d=NULL,f=NULL,fmod=NULL,c=NULL)
  {

### computes some forecast quality measures for some alternative forecasts
 
### requires "forecast", "stats" and "xts" packages

### y - a numeric or a column matrix of a dependent variable,

### x - a matrix of independent variables (drivers), different columns correspond to different independent variables

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

### c - a parameter indicating whether constant is included



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
if (is.null(f)) { f <- rep(TRUE,10) }
if (is.null(c)) { c <- TRUE }

y <- as.matrix(y)
x <- as.matrix(x)

######################### naive
##################################################

if (f[1]==TRUE)
{

y.naive <- (c(NA,as.vector(y)))[-(1+length(as.vector(y)))]
coeff.naive <- matrix(NA,nrow=1,ncol=ncol(x)+1)
pval.naive <- matrix(NA,nrow=1,ncol=ncol(x)+1)

}
else
{
y.naive <- NULL
coeff.naive <- NULL
pval.naive <- NULL
}

######################### OLS
##################################################

if (f[2]==TRUE)
{

if (c==TRUE)
  {
    m <- lm(y[1:(initial.period-1)] ~ x[1:(initial.period-1),,drop=FALSE])
    xx <- cbind(1,x)
  }
else
  {
    m <- lm(y[1:(initial.period-1)] ~ x[1:(initial.period-1),,drop=FALSE] -1)
    xx <- x
  }

y.ols <- as.vector(t(as.matrix(m$coefficients)) %*% t(xx))
mm <- summary(m)
if (all(is.finite(m$coefficients)) && all(is.finite(mm$coefficients[,4])))
  {
    coeff.ols <- t(mm$coefficients[,1,drop=FALSE])
    pval.ols <- t(mm$coefficients[,4,drop=FALSE])
    
    if (c==FALSE)
      {
        coeff.ols <- cbind(0,coeff.ols)
        pval.ols <- cbind(1,pval.ols)
      }
  }
else
  {
    coeff.ols <- matrix(NA,nrow=1,ncol=ncol(x)+1)
    pval.ols <- matrix(NA,nrow=1,ncol=ncol(x)+1)
  }
colnames(coeff.ols) <- c("const",colnames(x))
colnames(pval.ols) <- c("const",colnames(x))

}
else
{
y.ols <- NULL
coeff.ols <- NULL
pval.ols <- NULL
}

######################### recursive OLS
##################################################

if (f[3]==TRUE)
{

if (c==TRUE)
  {
    m <- rec.reg(y=y,x=x,c=TRUE)
  }
else
  {
    m <- rec.reg(y=y,x=x,c=FALSE)
  }
y.rec.ols <- as.vector(m$y.hat)
coeff.rec.ols <- m$coeff.
pval.rec.ols <- m$p.val.

if (c==FALSE)
  {
    coeff.rec.ols <- cbind(0,coeff.rec.ols)
    pval.rec.ols <- cbind(1,pval.rec.ols)
  } 
   
colnames(coeff.rec.ols) <- c("const",colnames(x))
colnames(pval.rec.ols) <- c("const",colnames(x))

}
else
{
y.rec.ols <- NULL
coeff.rec.ols <- NULL
pval.rec.ols <- NULL

}

######################### rolling OLS
##################################################

if (f[4]==TRUE)
{

if (c==TRUE)
  {
    m <- roll.reg(y=y,x=x,window=window,c=TRUE)
  }
else
  {
    m <- roll.reg(y=y,x=x,window=window,c=FALSE)
  }
y.roll.ols <- m$y.hat
coeff.roll.ols <- m$coeff.
pval.roll.ols <- m$p.val.

if (c==FALSE)
  {
    coeff.roll.ols <- cbind(0,coeff.roll.ols)
    pval.roll.ols <- cbind(1,pval.roll.ols)
  } 
    
colnames(coeff.roll.ols) <- c("const",colnames(x))
colnames(pval.roll.ols) <- c("const",colnames(x))

}
else
{
y.roll.ols <- NULL
coeff.roll.ols <- NULL
pval.roll.ols <- NULL

}

######################### TVP
##################################################

if (f[5]==TRUE)

{

if (c==TRUE)
  {
    m <- tvp(y=y,x=x,V=1,lambda=0.99,c=TRUE)
  }
else
  {
    m <- tvp(y=y,x=x,V=1,lambda=0.99,c=FALSE)
  }
y.tvp <- m$y.hat
coeff.tvp <- m$thetas
if (c==FALSE) { coeff.tvp <- cbind(0,coeff.tvp) }
pval.tvp <- matrix(NA,nrow=1,ncol=ncol(x)+1)
colnames(coeff.tvp) <- c("const",colnames(x))

}
else
{
y.tvp <- NULL
coeff.tvp <- NULL
pval.tvp <- NULL
}

######################### AR(1)
##################################################

if (f[6]==TRUE)
{

x1 <- (as.vector(y))[1:(length(as.vector(y))-1)]

yy <- (as.vector(y))[-1]

if (c==TRUE)
  {
    m <- lm(yy[1:(initial.period-1)] ~ x1[1:(initial.period-1)])
    xx <- cbind(1,x1)
  }
else
  {
    m <- lm(yy[1:(initial.period-1)] ~ x1[1:(initial.period-1)] -1)
    xx <- x1
  }
y.ar1 <- as.vector(t(as.matrix(m$coefficients)) %*% t(xx))
y.ar1 <- c(NA,y.ar1)

rm(x1,yy)

mm <- summary(m)

if (all(is.finite(m$coefficients)) && all(is.finite(mm$coefficients[,4])))
  {
    coeff.ar1 <- t(mm$coefficients[,1,drop=FALSE])
    pval.ar1 <- t(mm$coefficients[,4,drop=FALSE])
    
    if (c==FALSE)
      {
        coeff.ar1 <- cbind(0,coeff.ar1)
        pval.ar1 <- cbind(1,pval.ar1)
      }
  }
else
  {
    coeff.ar1 <- matrix(NA,ncol=2,nrow=1)
    pval.ar1 <- matrix(NA,ncol=2,nrow=1)
  }
colnames(coeff.ar1) <- c("const","ar1")
colnames(pval.ar1) <- c("const","ar1")
}
else
{
y.ar1 <- NULL
coeff.ar1 <- NULL
pval.ar1 <- NULL

}

######################### AR(2)
##################################################

if (f[7]==TRUE)
{

x1 <- (as.vector(y))[2:(length(as.vector(y))-1)]
x2 <- (as.vector(y))[1:(length(as.vector(y))-2)]
x2 <- cbind(x1,x2)

yy <- (as.vector(y))[-c(1,2)]

if (c==TRUE)
  {
    m <- lm(yy[1:(initial.period-1)] ~ x2[1:(initial.period-1),,drop=FALSE])
    xx <- cbind(1,x2)
  }
else
  {
    m <- lm(yy[1:(initial.period-1)] ~ x2[1:(initial.period-1),,drop=FALSE] -1)
    xx <- x2
  }
mm <- summary(m)

y.ar2 <- as.vector(t(as.matrix(m$coefficients)) %*% t(xx))
y.ar2 <- c(NA,NA,y.ar2)

rm(x1,yy,x2)

if (all(is.finite(m$coefficients)) && all(is.finite(mm$coefficients[,4])))
  {
    coeff.ar2 <- t(mm$coefficients[,1,drop=FALSE])
    pval.ar2 <- t(mm$coefficients[,4,drop=FALSE])
    
    if (c==FALSE)
      {
        coeff.ar2 <- cbind(0,coeff.ar2)
        pval.ar2 <- cbind(1,pval.ar2)
      }
  }
else
  {
    coeff.ar2 <- matrix(NA,ncol=3,nrow=1)
    pval.ar2 <- matrix(NA,ncol=3,nrow=1)
  }
colnames(coeff.ar2) <- c("const","ar1","ar2")
colnames(pval.ar2) <- c("const","ar1","ar2")

}
else
{
y.ar2 <- NULL
coeff.ar2 <- NULL
pval.ar2 <- NULL
}

######################### auto ARIMA
##################################################

if (f[8]==TRUE)
{

if (c==TRUE)
  {
    m <- auto.arima(as.vector(y)[1:(initial.period-1)],allowmean=TRUE)
  }
else
  {
    m <- auto.arima(as.vector(y)[1:(initial.period-1)],allowmean=FALSE)
  }
y.auto.arima <- y[1:(initial.period-1)] - m$residuals
if (initial.period==1)
  {
    y.auto.arima <- c(as.vector(y.auto.arima),as.vector(forecast(object=m,h=length(y)-1)$mean)) 
  }
else
  {
    y.auto.arima <- c(as.vector(y.auto.arima),as.vector(forecast(object=m,h=length(y)-(initial.period-1))$mean)) 
  }
  
se <- (diag(m$var.coef))^0.5

coeff.auto.arima <- m$coef
if (any(dim(se)==0))
  {
    pval.auto.arima <- 1
  }
else
  {
    pval.auto.arima <- (1-pnorm(abs(coeff.auto.arima)/se))*2
  }

if (c("intercept") %in% names(coeff.auto.arima))
  {
    j <- which(names(coeff.auto.arima)=="intercept")
    names(coeff.auto.arima)[j] <- "const"
  }
coeff.auto.arima <- t(as.matrix(coeff.auto.arima))
pval.auto.arima <- t(as.matrix(pval.auto.arima))
colnames(pval.auto.arima) <- colnames(coeff.auto.arima)
if (c("const") %in% colnames(coeff.auto.arima))
  {
    j <- which(colnames(coeff.auto.arima)=="const")
    coeff.auto.arima <- cbind(coeff.auto.arima[,j],coeff.auto.arima[,-j])
    pval.auto.arima <- cbind(pval.auto.arima[,j],pval.auto.arima[,-j])
  }

}
else
{
y.auto.arima <- NULL
coeff.auto.arima <- NULL
pval.auto.arima <- NULL
}

######################### TVP-AR(1)
##################################################

if (f[9]==TRUE)
{

x1 <- (as.vector(y))[1:(length(as.vector(y))-1)]
yy <- (as.vector(y))[-1]

x1 <- as.matrix(x1)
yy <- as.matrix(yy)
colnames(x1) <- c("ar1")

if (c==TRUE)
  {
    m <- tvp(y=yy,x=x1,V=1,lambda=0.99,c=TRUE)
  }
else
  {
    m <- tvp(y=yy,x=x1,V=1,lambda=0.99,c=FALSE)
  }
y.tvp.ar1 <- m$y.hat
y.tvp.ar1 <- c(NA,y.tvp.ar1)

rm(x1,yy)

coeff.tvp.ar1 <- m$thetas
if (c==FALSE) { coeff.tvp.ar1 <- cbind(0,coeff.tvp.ar1) }
coeff.tvp.ar1 <- rbind(rep(NA,2),coeff.tvp.ar1)
pval.tvp.ar1 <- matrix(NA,nrow=1,ncol=2)
colnames(coeff.tvp.ar1) <- c("const","ar1")

}
else
{
y.tvp.ar1 <- NULL
coeff.tvp.ar1 <- NULL
pval.tvp.ar1 <- NULL

}

######################### TVP-AR(2)
##################################################

if (f[10]==TRUE)
{

x1 <- (as.vector(y))[2:(length(as.vector(y))-1)]
x2 <- (as.vector(y))[1:(length(as.vector(y))-2)]
x2 <- cbind(x1,x2)

yy <- (as.vector(y))[-c(1,2)]

x2 <- as.matrix(x2)
yy <- as.matrix(yy)
colnames(x2) <- c("ar1","ar2")

if (c==TRUE)
  {
    m <- tvp(y=yy,x=x2,V=1,lambda=0.99,c=TRUE)
  }
else
  {
    m <- tvp(y=yy,x=x2,V=1,lambda=0.99,c=FALSE)
  }

y.tvp.ar2 <- m$y.hat
y.tvp.ar2 <- c(NA,NA,y.tvp.ar2)

rm(x1,yy,x2)

coeff.tvp.ar2 <- m$thetas
if (c==FALSE) { coeff.tvp.ar2 <- cbind(0,coeff.tvp.ar2) }
coeff.tvp.ar2 <- rbind(rep(NA,3),rep(NA,3),coeff.tvp.ar2)

pval.tvp.ar2 <- matrix(NA,nrow=1,ncol=3)
colnames(coeff.tvp.ar2) <- c("const","ar1","ar2")

}
else
{
y.tvp.ar2 <- NULL
coeff.tvp.ar2 <- NULL
pval.tvp.ar2 <- NULL

}

##################################################

fq <- list(y.naive,
           y.ols,
           y.rec.ols,
           y.roll.ols,
           y.tvp,
           y.ar1,
           y.ar2,
           y.auto.arima,
           y.tvp.ar1,
           y.tvp.ar2)

coeff <- list(coeff.naive,
              coeff.ols,
              coeff.rec.ols,
              coeff.roll.ols,
              coeff.tvp,
              coeff.ar1,
              coeff.ar2,
              coeff.auto.arima,
              coeff.tvp.ar1,
              coeff.tvp.ar2)

pval <- list(pval.naive,
             pval.ols,
             pval.rec.ols,
             pval.roll.ols,
             pval.tvp,
             pval.ar1,
             pval.ar2,
             pval.auto.arima,
             pval.tvp.ar1,
             pval.tvp.ar2)

          
fq2 <- fq[f]

for (i in 1:10)
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

rnames <- c("naive","OLS","rec. OLS","roll. OLS","TVP","AR(1)","AR(2)","auto ARIMA","TVP-AR(1)","TVP-AR(2)")
rownames(fq) <- rnames[f]
names(fq2) <- rnames[f]
colnames(fq) <- c("ME","RMSE","MAE","MPE","MAPE","HR")

coeff <- coeff[f]
pval <- pval[f]
names(coeff) <- rnames[f]
names(pval) <- rnames[f]


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

r <- list(round(fq,digits=4),fq2,as.matrix(y),coeff,pval)
names(r) <- c("summary","y.hat","y","coeff.","p.val.")
class(r) <- "altf"
return(r)

  }
  