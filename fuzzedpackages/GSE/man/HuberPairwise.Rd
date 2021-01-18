\name{HuberPairwise}
\alias{HuberPairwise}
\title{ Quadrant Covariance and Huberized Pairwise Scatter }
\description{
  Computes the Quadrant Covariance (QC) or Huberized Pairwise Scatter as described in 
  Alqallaf et al. (2002). 
}
\usage{ 
    HuberPairwise( x, psi=c("huber","sign"), c0=1.345, computePmd=TRUE)
}
\arguments{
  \item{x}{a matrix or data frame. May contain missing values, but cannot contain columns with completely missing entries. }
  \item{psi}{loss function to be used in computing pairwise scatter. Default is "huber". If \code{psi="sign"}, this yields QC. Other value includes "huber". }
  \item{c0}{tuning constant for the huber function. \code{c0=0} would yield QC. Default is \code{c0=1.345}. This parameter is unnecessary if \code{psi='sign'}.}
  \item{computePmd}{logical indicating whether to compute partial Mahalanobis distances (pmd) and adjusted pmd. Default is \code{TRUE}. }
}
\details{
   As described in Alqallaf et al. (2002), this estimator requires a robust scale estimate and a location M-estimate, 
   which will be used to transform the data through a loss-function to be outlier-free. Currently, this function takes
   MADN (normalized MAD) and median as the robust scale and location estimate to save computation time. By default, the 
   loss function \code{psi} is a sign function, but users are encouraged to also try Huberized scatter with the loss function 
   as \eqn{\psi_c(x) = min( max(-c, x), c), c > 0, c=1.345}. The function does not adjust for intrinsic bias as described in 
   Alqallaf et al. (2002). Missing values will be replaced by the corresponding column's median.
}
\value{
  An S4 object of class \code{\link{HuberPairwise-class}} which is a subclass of the virtual class \code{\link{CovRobMiss-class}}.
  The output S4 object contains the following slots:

  \tabular{ll}{
    \code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
    \code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
    \code{pmd} \tab Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. \cr
    \code{pmd.adj} \tab Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. \cr
    \code{pu} \tab Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. \cr
    \code{R} \tab Estimated correlation matrix. Not meant to be accessed. \cr
    \code{call} \tab Object of class \code{"language"}. Not meant to be accessed. \cr
    \code{x} \tab Input data matrix. Not meant to be accessed. \cr
    \code{p} \tab Column dimension of input data matrix. Not meant to be accessed. \cr
    \code{estimator}\tab Character string of the name of the estimator used. Not meant to be accessed. \cr
    }
}  
\references{
  Alqallaf, F.A., Konis, K. P., R. Martin, D., Zamar, R. H. (2002).
  Scalable Robust Covariance and Correlation Estimates for Data Mining. In Proceedings 
  of the Seventh ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. Edmonton.
}
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }


