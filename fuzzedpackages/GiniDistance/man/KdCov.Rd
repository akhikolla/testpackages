\name{KdCov}
\title{Kernel Distance Covariance Statistics}
\alias{KdCov}
\description{
  Computes Kernel distance covariance statistics,
  in which Xs are quantitative, Y are categorical, sigma is kernel standard deviation and returns the measures of dependence.
}
\usage{
  KdCov(x, y, sigma)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{sigma}{ kernel standard deviation}
}
\details{
  \code{KdCov} compute distance correlation statistics.
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels. 
  
  Distance covariance was introduced in (\cite{Szekely07}) as a dependence measure between random variables \eqn{X \in {R}^p} and \eqn{Y \in {R}^q}. If \eqn{X} and \eqn{Y} are embedded into RKHS's induced by \eqn{\kappa_X} and \eqn{\kappa_Y}, respectively, the generalized distance covariance of \eqn{X} and \eqn{Y} is (\cite{Sejdinovic13}):
  \deqn{\begin{array}{c}
    \mathrm{dCov}_{\kappa_X,\kappa_Y}(X,Y) = {E}d_{\kappa_X}(X,X^{\prime})d_{\kappa_Y}(Y,Y^{\prime}) + {E}d_{\kappa_X}(X,X^{\prime}){E}d_{\kappa_Y}(Y,Y^{\prime}) \\
    - 2{E}\left[{E}_{X^{\prime}}d_{\kappa_X}(X,X^{\prime}) {E}_{Y^{\prime}}d_{\kappa_Y}(Y,Y^{\prime})\right].\label{dCovkk}
    \end{array}
  }
  In the case of \eqn{Y} being categorical, one may embed it using a set difference kernel \eqn{\kappa_Y},
  \deqn{\label{setdiff}
    \kappa_Y(y,y^{\prime}) = \left\{ \begin{array}{cc}
    \frac{1}{2} &  if \;y = y^{\prime},\\ 0 & otherwise.
    \end{array} \right.
  }
  This is equivalent to embedding \eqn{Y} as a simplex with edges of unit length (\cite{Lyons13}), i.e., \eqn{L_k} is represented by a \eqn{K} dimensional vector of all zeros except its \eqn{k}-th dimension, which has the value \eqn{\frac{\sqrt{2}}{2}}. 
  The distance induced by \eqn{\kappa_Y} is called the set distance, i.e., \eqn{d_{\kappa_Y}(y,y^{\prime})=0} if \eqn{y=y^{\prime}} and \eqn{1} otherwise. Using the set distance, we have the following results on the generalized distance covariance between a numerical 
  and a categorical random variable.
  \deqn{\mathrm{dCov}_{\kappa_X,\kappa_Y}(X,Y) := \mathrm{dCov}_{\kappa_X}(X,Y) \nonumber = \sum_{k=1}^{K} p_k^2 \left[2 {E}d_{\kappa_X}(X_k,X) - {E}d_{\kappa_X}(X_k,{X_k}^{\prime}) - {E}d_{\kappa_X}(X,X^{\prime}) \right].\label{dCovk}}
    
}           
  
\references{
  
  Sejdinovic, D., Sriperumbudur, B., Gretton, A.  and Fukumizu, K. (2013). Equivalence of Distance-based and RKHS-based Statistics in Hypothesis Testing, The Annals of Statistics, 41 (5),  2263-2291. 
    
  Zhang, S., Dang, X., Nguyen, D. and Chen, Y. (2019). Estimating feature - label dependence using Gini distance statistics. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence (submitted)}.
    
}
\value{
  \code{KdCov} returns the sample kernel distance correlation
}
\seealso{
  \code{\link{KgCov}}  \code{\link{KgCor}}  \code{\link{dCov}}
}
\examples{
  x<-iris[,1:4]
  y<-unclass(iris[,5])
  KdCov(x, y, sigma=1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
