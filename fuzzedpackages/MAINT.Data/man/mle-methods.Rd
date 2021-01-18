\name{mle-methods}
\docType{methods}
\alias{mle-methods}
\alias{mle,IData-method}
\alias{mle}   
\title{ Methods for function mle in Package `MAINT.Data'}

\description{Performs maximum likelihood estimation for parametric models of interval data}

\usage{

\S4method{mle}{IData}(Idt, Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4,
  SelCrit=c("BIC","AIC"), k2max=1e6, OptCntrl=list(), \dots)

}

\arguments{
  \item{Idt}{An IData object representing interval-valued entities.}

  \item{Model}{The joint distribution assumed for the MidPoint and LogRanges. Current alternatives are \dQuote{Normal} for Gaussian,
distributions, \dQuote{SNNormal} for Skew-Normal and \dQuote{NrmandSKN} for both Gaussian and Skew-Normal distributions.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

   \item{SelCrit}{The model selection criterion.}

   \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

   \item{OptCntrl}{List of optional control parameters to passed to the optimization routine. See the documentation of RepLOptim for a description of the available options.}

  \item{\dots}{Other named arguments.}
}

\references{
Azzalini, A. and Dalla Valle, A. (1996), The multivariate skew-normal distribution. \emph{Biometrika} \bold{83}(4), 715--726.\cr

Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.
}

\seealso{
\code{\link{IData}}, \code{\link{RepLOptim}}
}

\keyword{methods}
\keyword{maximum likelihood estimation}
\keyword{interval data}

\examples{
# Create an Interval-Data object containing the intervals of temperatures by quarter 
# for 899 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8])

# Estimate parameters by maximum likelihood assuming a Gaussian distribution

ChinaE <- mle(ChinaT)
cat("China maximum likelhiood estimation results =\n")
print(ChinaE)
cat("Standard Errors of Estimators:\n")
print(stdEr(ChinaE))
}
