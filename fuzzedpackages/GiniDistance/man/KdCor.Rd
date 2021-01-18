\name{KdCor}
\title{Kernel Distance Correlation Statistics}
\alias{KdCor}
\description{
  Computes Kernel distance correlation statistics,
  in which Xs are quantitative, Y are categorical, sigma is kernel standard deviation and returns the measures of dependence.
}
\usage{
  KdCor(x, y, sigma)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{sigma}{ kernel standard deviation}
}
\details{
  \code{KdCor} compute distance correlation statistics.
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels. 
  
  The kernel distance correlation is defined as follow.
  
  \deqn{{dCor}_{\kappa_X,\kappa_Y}(\mathbf{X}, Y) = \frac{\mbox{ dCov}_{\kappa_X,\kappa_Y}(\mathbf{X}, Y)}{ \sqrt{\mbox{ dCov}_{\kappa_X,\kappa_X}(\mathbf{X},\mathbf{X})} \sqrt{\mbox{ dCov}_{\kappa_Y,\kappa_Y}(Y,Y)}}}
  
  where
  
  \deqn{\begin{array}{c}
  \mathrm{dCov}_{\kappa_X,\kappa_Y}(X,Y) = {E}d_{\kappa_X}(X,X')d_{\kappa_Y}(Y,Y') + {E}d_{\kappa_X}(X,X'){E}d_{\kappa_Y}(Y,Y') \\
  - 2{E}\left[{E}_{X'}d_{\kappa_X}(X,X') {E}_{Y'}d_{\kappa_Y}(Y,Y')\right].\label{dCovkk}
  \end{array}
 }
    
}
\value{
\code{KdCor} returns the sample kernel distance correlation
}
    
\references{

  Sejdinovic, D., Sriperumbudur, B., Gretton, A.  and Fukumizu, K. (2013). Equivalence of Distance-based and RKHS-based Statistics in Hypothesis Testing, The Annals of Statistics, 41 (5),  2263-2291. 
  
  Zhang, S., Dang, X., Nguyen, D. and Chen, Y. (2019). Estimating feature - label dependence using Gini distance statistics. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence (submitted)}.

      
}
\seealso{
  \code{\link{KgCov}}  \code{\link{KgCor}}  \code{\link{dCor}}
}
\examples{
  x<-iris[,1:4]
  y<-unclass(iris[,5])
  KdCor(x, y, sigma=1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
