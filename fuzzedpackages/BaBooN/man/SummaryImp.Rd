\name{summary.imp}

\alias{summary.imp}

\title{Summary method for objects of class \sQuote{imp}}

\description{Returns some information about the incomplete data set and the imputation process.} 

\usage{\method{summary}{imp}(object,...)}

\arguments{
  \item{object}{Either with \code{BBPMM} or \code{BBPMM.row} generated object.}
  \item{...}{Arguments to be passed to or from other functions.}
}

\details{
Returns information about the percentage of missing data as well as about the imputation variant, 
the number of (multiple) imputations and the number of iterations between two imputations.
}

\author{Florian Meinfelder}

\seealso{\code{\link{BBPMM}}, \code{\link{BBPMM.row}}}

\examples{
### sample data set with non-normal variables and two different
### missingness patterns
n <- 50
x1 <- round(runif(n,0.5,3.5))
x2 <- as.factor(c(rep(1,10),rep(2,25),rep(3,15)))
x3 <- round(rnorm(n,0,3))
y1 <- round(x1-0.25*(x2==2)+0.5*x3+rnorm(n,0,1))
y1 <- ifelse(y1<1,1,y1)
y1 <- as.factor(ifelse(y1>4,5,y1))
y2 <- x1+rnorm(n,0,0.5)
y3 <- round(x3+rnorm(n,0,2))
data1 <- as.data.frame(cbind(x1,x2,x3,y1,y2,y3))
misrow1 <- sample(n,20)
misrow2 <- sample(n,15)
misrow3 <- sample(n,10)
is.na(data1[misrow1, 4]) <- TRUE
is.na(data1[misrow2, 5]) <- TRUE
is.na(data1[misrow2, 6]) <- TRUE

### imputation
imputed.data <- BBPMM(data1, nIter=5, M=5)
summary(imputed.data)
}
