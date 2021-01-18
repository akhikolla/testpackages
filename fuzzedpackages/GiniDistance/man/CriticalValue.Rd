\name{CriticalValue}
\title{Find a critical value by permutation test of dependence between X and Y using kernel (Gini) distance covariance or correlation statistics}
\alias{CriticalValue}
\description{
  Find a critical value by permutation test using variance of kernel (Gini) distance covariance or correlation statistics,
  in which Xs are quantitative, Y are categorical, sigma is kernel standard deviation, alpha is an exponent on Euclidean distance and returns the critical value of the measures of dependence.
}
\usage{
  CriticalValue(x, y, sigma, alpha, level, M = 1000, method)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{sigma}{ kernel standard deviation}
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
  \item{level}{ significance level of the test, the default value = 0.05}
  \item{M}{ number of permutations}
  \item{method}{ string name of the method for permutation test, e.g. gCov}
}
\details{
  \code{CriticalValue} compute the critical value of a dependence test of a kernel (Gini) distance covariance or correlation statistics.
  It is a self-contained R function returning the critical value of the measure of dependence statistics. 
  
  The critical value of the test of significance level \eqn{\gamma}, however, is obtained by a permutation procedure.
  Let \eqn{\nu = 1: n} be the vector of original sample indices of the sample for \eqn{Y} labels and \eqn{\hat{\rho}_g(\alpha) = \hat{\rho}(\nu;\alpha)}. 
  Let \eqn{\pi(\nu)} denote a permutation of the elements of \eqn{\nu} and the corresponding \eqn{\hat{\rho}_g(\pi;\alpha)} is computed. 
  Under the \eqn{{\cal H}_0}, \eqn{\hat{\rho}_g(\nu)}
  and \eqn{\hat{\rho}_g(\pi;\alpha)} are identically distributed for every permutation \eqn{\pi} of \eqn{\nu}. 
  Hence, based on \eqn{M} permutations,  the critical value \eqn{q_{\gamma}} is estimated by the \eqn{(1-\gamma)100\%} sample 
  quantile of \eqn{\hat{\rho}_g(\pi_m;\alpha)}, \eqn{m=1,...,M}. Usually \eqn{100\leq M\leq 1000} is sufficient 
  for a good estimation on the critical value.   
    
  See \code{\link{PermutationTest}} for a test of multivariate independence
  based on the (Gini) distance statistic.
}
\value{
  \code{CriticalValue} returns return the critical value of the measures of the dependence of the permutation test of a specified function
}
\seealso{
  \code{\link{PermutationTest}} 
}
\examples{
  n = 50
  x <- runif(n)
  y <- c(rep(1,n/2),rep(2,n/2))
  CriticalValue(x, y, sigma=1, alpha=2, level=0.04, M = 1000, method='KgCov') 
}
\keyword{ multivariate }
\concept{ independence }
\concept{ permutation test}
\concept{ Gini statistics }
