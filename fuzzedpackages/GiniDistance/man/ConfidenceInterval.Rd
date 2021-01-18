\name{ConfidenceInterval}
\title{Confidence Interval of Dependence measure}
\alias{ConfidenceInterval}
\description{
  Find confidence intervals for dependence measures in which Xs are quantitative, Y are categorical using jack-knife method.
}
\usage{
  ConfidenceInterval(x, y, sigma, alpha, level, method)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{sigma}{ kernel parameter}
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
  \item{level}{ level of confidence, in [0,1]}
  \item{method}{ name of dependence measure which can chosen from "gCor","gCov","dCor","dCov","KgCor", "KgCov", "KdCor" and "KdCov"}
}
\details{
  \code{ConfidenceInterval} compute the confidence interval of the distance correlation statistics.
  It is a self-contained R function returning a variance of the measure of dependence statistics. 
  
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels. \code{alpha} if missing by default is 1, otherwise it is exponent on the Euclidean distance.
  
  Suppose a sample data \eqn{ {\mathcal{D}} =\{(\boldmath{x}_i,y_i)\} } for \eqn{i = 1,...,n} available. The confidence interval is built upon the asymptotic normality of sample dependence statistic. The asymptotic variance is estimated by the Jackknife method. 
  More details refer to Shao and Tu (1996).   
}
\value{
  \code{ConfidenceInterval} returns the confidence interval of distance correlation
}

\references{
  
  Dang, X., Nguyen, D., Chen, Y. and Zhang, J. (2019). A new Gini correlation between quantitative and qualitative variables. Submitted. 
  
  Shao, J. and Tu, D. (1996). {The Jackknife and Bootstrap}. Springer, New York.
} 

\examples{
  x <- iris[,1:4]
  y <- unclass(iris[,5])
  ConfidenceInterval(x, y, alpha=1, level=0.95, method='gCor')
}

\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
