\name{MI.inference}

\alias{MI.inference}

\title{Multiple Imputation inference}

\description{\sQuote{MI.inference} applies Rubin's combining rules to
  estimated quantities of interest that are based on multiply imputed
  data sets. The function requires as input two vectors of length M for the
  estimate and its variance.  
}

\usage{MI.inference(thetahat, varhat.thetahat, alpha=0.05)}

\arguments{
  \item{thetahat}{A vector of length M containing estimates of the
    quantity of interest based on multiply imputed data sets.}
  \item{varhat.thetahat}{A vector of length M containing the
    corresponding variances of \code{thetahat}.}
  \item{alpha}{The significance level at which lower and upper bound are
    calculated. DEFAULT=0.05}
}

\value{
  \item{MI.Est}{A scalar containing the MI estimate of the quantity of
    interest (i.e. an estimator averaged over all M data sets).}
  \item{MI.Var}{The Multiple Imputation variance.}
  \item{CI.low}{The lower bound of the MI confidence interval.}
  \item{CI.up}{The upper bound of the MI confidence interval.}
  \item{BVar}{The estimated \emph{between variance}.}
  \item{WVar}{The estimated \emph{within variance}.}
}

\details{Multiple Imputation (Rubin, 1987) of missing data is a
  generally accepted way to get correct variance estimates for a
  particular quantity of interest in the presence of missing
  data. \code{MI.inference} estimates the \emph{within variance} W and
  \emph{between variance} B, and combines them to the  \emph{total
    variance} T. Based on the output, further analysis figures, such as
  the \emph{fraction of missing information} can be calculated.}


\references{
  Rubin, D.B. (1987) \emph{Multiple Imputation for Non-Response
  in Surveys}. New York: John Wiley & Sons, Inc.
}

\examples{
\dontrun{
### example 1
n <- 100
x1 <- round(runif(n,0.5,3.5))
x2 <- round(runif(n,0.5,4.5))
x3 <- runif(n,1,6)
y1 <- round(x1-0.25*x2+0.5*x3+rnorm(n,0,1))
y1 <- ifelse(y1<2,2,y1)
y1 <- as.factor(ifelse(y1>4,5,y1))
y2 <- x3+rnorm(n,0,2)
y3 <- as.factor(ifelse(x2+rnorm(n,0,2)>2,1,0))
mis1 <- sample(100,20)
mis2 <- sample(100,30)
mis3 <- sample(100,25)
data1 <- data.frame("x1"=x1,"x2"=x2,"x3"=x3,
                    "y1"=y1,"y2"=y2,"y3"=y3)
is.na(data1$y1[mis1]) <- TRUE
is.na(data1$y2[mis2]) <- TRUE
is.na(data1$y3[mis3]) <- TRUE
imputed.data <- BBPMM(data1, M=5, nIter=5)

MI.m.meany2.hat <- sapply(imputed.data$impdata,
                          FUN=function(x) mean(x$y2))
                          
MI.v.meany2.hat <- sapply(imputed.data$impdata,
                          FUN=function(x) var(x$y2)/length(x$y2))

### MI inference
MI.y2 <- MI.inference(MI.m.meany2.hat,
                      MI.v.meany2.hat, alpha=0.05)

MI.y2$MI.Est
MI.y2$MI.Var


################################################################
### example 2: a small simulation example

### simple additional function to calculate coverages:         #

coverage <- function(value, bounds) {
  ifelse(min(bounds) <= value && max(bounds) >= value, 1, 0)
}
### value            : true value                              #
### bounds           : vector with two elements (upper and     #
###                    lower bound of the CI)                  #

### sample size
n <- 100
### true value for the mean of y2
m.y2 <- 3.5
y2.cover <- vector(length=n)
set.seed(1000)

### 100 data generations
time1 <- Sys.time()
for (i in 1:100) {
  x1 <- round(runif(n,0.5,3.5))
  x2 <- round(runif(n,0.5,4.5))
  x3 <- runif(n,1,6)
  y1 <- round(x1-0.25*x2+0.5*x3+rnorm(n,0,1))
  y1 <- ifelse(y1<2,2,y1)
  y1 <- as.factor(ifelse(y1>4,5,y1))
  y2 <- x3+rnorm(n,0,2)
  y3 <- as.factor(ifelse(x2+rnorm(n,0,2)>2,1,0))
  mis1 <- sample(n,20)
  mis2 <- sample(n,30)
  mis3 <- sample(n,25)
  data1 <- data.frame("x1"=x1,"x2"=x2,"x3"=x3,
                      "y1"=y1,"y2"=y2,"y3"=y3)
  is.na(data1$y1[mis1]) <- TRUE
  is.na(data1$y2[mis2]) <- TRUE
  is.na(data1$y3[mis3]) <- TRUE
  
  sim.imp <- BBPMM(data1, M=3, nIter=2,
                   stepmod="", verbose=FALSE)

  MI.m.meany2.hat <- sapply(sim.imp$impdata,
                            FUN=function(x) mean(x$y2))

  MI.v.meany2.hat <- sapply(sim.imp$impdata,
                            FUN=function(x) 
                            var(x$y2)/length(x$y2))
### MI inference
  MI.y2 <- MI.inference(MI.m.meany2.hat, MI.v.meany2.hat,
                        alpha=0.05)

  y2.cover[i] <- coverage(m.y2, c(MI.y2$CI.low,MI.y2$CI.up))
}
time2 <- Sys.time()
difftime(time2, time1, unit="secs")

### coverage estimator (alpha=0.05):
mean(y2.cover)

}
}
