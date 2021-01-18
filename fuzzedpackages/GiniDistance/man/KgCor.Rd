\name{KgCor}
\title{ Kernel Gini Distance Correlation Statistics}
\alias{KgCor}
\description{
  Computes Kernel Gini distance correlation statistics,
  in which Xs are quantitative, Y are categorical, sigma is kernel standard deviation, alpha is an exponent on the Euclidean distance and returns the kernel Gini mean difference.
}
\usage{
  KgCor(x, y, sigma)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{sigma}{ kernel standard deviation}
}
\details{
  \code{Kgcor} compute kernel Gini distance correlation statistics for data.
  It is a self-contained R function dealing with both univariate and multivariate data. 
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels.
  
  Gini distance correlation are generalized to RKHS, \eqn{\mathcal{H}_\kappa}, as
    
  \deqn{ 
  \mathrm{gCor}_\kappa(X,Y)  = \frac{\sum_{k=1}^{K} p_k \left[ 2 {E}d_\kappa(X_k,X) - 
  {E}d_\kappa(X_k,{X_k}') - {E}d_\kappa(X,X')\right]}{{E}d_\kappa(X,X')}.\label{gCork}}
  In this case, we use the default Gaussian distance function
  \deqn{
  d_\kappa(x,x') = \sqrt{1-e^{-\frac{|x-x'|_q^2}{\sigma^2}}},
  }
  induced by a weighted Gaussian kernel, \eqn{\kappa(x,x') = \frac{1}{2}e^{-\frac{|x-x'|_q^2}{\sigma^2}}.} 
        
}
\value{
  \code{KgCor} returns the sample Kernel Gini distance correlation between \code{x} and \code{y}.
}
\seealso{
  \code{\link{gCov}}  \code{\link{gCor}}  \code{\link{dCor}}
}
        
\references{
  Zhang, S., Dang, X., Nguyen, D. and Chen, Y. (2019). Estimating feature - label dependence using Gini distance statistics. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence (submitted)},
  \url{https://arXiv.org/pdf/1906.02171.pdf}
            
} 
\examples{
  x<-iris[,1:4]
  y<-unclass(iris[,5])
  KgCor(x, y, sigma=1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
