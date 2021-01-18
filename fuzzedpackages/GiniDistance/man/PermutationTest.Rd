\name{PermutationTest}
\title{Permutation test of dependence between X and Y using (Gini) distance covariance or correlation statistics}
\alias{PermutationTest}
\description{
  Perform permutation test using various dependence measures, 
  in which Xs are quantitative, Y are categorical, alpha is an exponent on Euclidean distance, sigma is kernel parameter in kernel methods and return the test statistic, critical value, p-value and decision of the test.
}
\usage{
  PermutationTest(x, y,  method, sigma, alpha, M = 200, level = 0.05)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{method}{ name of permutation test method and is chosen from one of the method list: dCov, dCor, KdCov, KdCor, gCov, gCor, KgCov, Kgcor}
  \item{sigma}{ kernel parameter for kenerl methods}
  \item{alpha}{ exponent on Euclidean distance, in (0,2), the default value = 1}
  \item{M}{number of permutations}  
  \item{level}{ significance level of the test, the default value = 0.05}
  
}
\details{
  \eqn{H_0:} X and Y are independent \eqn{\Longleftrightarrow H_0: F(x|y=1)=F(x|Y=2)=...=F(x|Y=K)}
  
  \code{PermutationTest} compute the p-value value of a permutation test of a (Gini) distance covariance or correlation statistics.
  It is a self-contained R function  the measure of dependence statistics. 
  
  The p-value is obtained by a permutation procedure.
  Let \eqn{\hat{\rho}(\nu)} be the sample dependnce measure based on the orginal sample indexed by \eqn{\nu=\{1,2,...,n\}}.  Let \eqn{\pi(\nu)} denote a permutation of the elements of \eqn{\nu} and the corresponding \eqn{\hat{\rho}(\pi)} is computed for the permutated data on y labels.   
  Under the \eqn{{\cal H}_0}, \eqn{\hat{\rho}(\nu)}
  and \eqn{\hat{\rho}(\pi)} are identically distributed for every permutation \eqn{\pi} of \eqn{\nu}. 
  Hence, based on \eqn{M} permutations,  the critical value \eqn{q_{\gamma}} is estimated by the \eqn{(1-\gamma)100\%} sample 
  quantile of \eqn{\hat{\rho}(\pi_m)}, \eqn{m=1,...,M} and the p-value is estimated by the proportion of
  \eqn{\hat{\rho}(\pi_m)} greater than \eqn{\hat{\rho}(\nu)}. Usually \eqn{100\leq M\leq 1000} is sufficient for a good estimation on the critical value or p-value.  The default value is \eqn{M=200}.
}    
    
\value{
  \code{PermutationTest} returns the p-value,  critical value and decision of the permutation test of a specified method.
}
    
\seealso{
  \code{\link{gCor}} \code{\link{gCov}} \code{\link{dCor}} \code{\link{dCov}} \code{\link{KgCov}} \code{\link{KgCov}} \code{\link{KdCov}} 
}
\examples{
  n = 50
  x <- runif(n)
  y <- c(rep(1,n/2),rep(2,n/2))
  PermutationTest(x, y, method = "gCor", alpha = 2, M = 50 )
}
\keyword{ multivariate }
\concept{ independence }
\concept{ permutation test}
\concept{ Gini statistics }
