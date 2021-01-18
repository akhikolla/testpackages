\name{rec.reg}
\alias{rec.reg}

\title{Computes Recursive Regression.}

\description{
This function computes Recursive Regression.  
}

\details{It might happen during computations that \code{\link[stats]{lm}} (which is used inside \code{rec.reg}) will produce \code{\link[base]{NA}} or \code{\link[base]{NaN}}. In such a case regression coefficients for a given period are taken as \eqn{0} and p-values for t-test for statistical significance of regression coefficients are taken as \eqn{1}.

It is not possible to set \code{c=FALSE} if \code{x=NULL}. In such a case the function will automatically reset \code{c=TRUE} inside the code.}

\usage{
rec.reg(y,x=NULL,c=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link[base]{matrix}} of a dependent variable}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns should correspond to different variables, if not specified only a constant will be used}
\item{c}{optional, \code{\link[base]{logical}}, a parameter indicating whether constant is included, if not specified \code{c=TRUE} is used, i.e., constant is included}
}

\value{class \code{reg} object, \code{\link[base]{list}} of
\item{$y.hat}{fitted (forecasted) values}
\item{$AIC}{Akaike Information Criterion (from the current set of observations)}
\item{$AICc}{Akaike Information Criterion with a correction for finite sample sizes (from the current set of observations)}
\item{$BIC}{Bayesian Information Criterion (from the current set of observations)} 
\item{$MSE}{Mean Squared Error (from the current set of observations)}
\item{$coeff.}{regression coefficients}
\item{$p.val}{p-values for t-test for statistical significance of regression coefficients}
\item{$y}{\code{y}, forecasted time-series}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
rec1 <- rec.reg(y=ld.wti,x=ld.drivers)
rec2 <- rec.reg(y=ld.wti)
}

\seealso{\code{\link{print.reg}}, \code{\link{summary.reg}}, \code{\link{plot.reg}}.}
