\name{emve}
\alias{emve}
\title{Extended Minimum Volume Ellipsoid (EMVE) in the presence of missing data}
\description{
  Computes the Extended S-Estimate (ESE) version of the minimum volume ellipsoid (EMVE), which is
  used as an initial estimator in Generlized S-Estimator (GSE) for missing data by default. 
}
\usage{emve(x, maxits=5, sampling=c("uniform","cluster"), n.resample, n.sub.size, seed) }
\arguments{
  \item{x}{a matrix or data frame. May contain missing values, but cannot contain columns with completely missing entries.}
  \item{maxits}{integer indicating the maximum number of iterations of Gaussian MLE calculation for each subsample. Default is 5.}
  \item{sampling}{which sampling scheme is to use: 'uniform' or 'cluster' (see Leung and Zamar, 2016). Default is 'uniform'.}
  \item{n.resample}{integer indicating the number of subsamples. Default is 15 for clustering-based subsampling and 500 for uniform subsampling.}
  \item{n.sub.size}{integer indicating the sizes of each subsample. Default is 2(p+1)/a for clustering-based subsampling and (p+1)/a for uniform subsampling, where a is proportion of non-missing cells.}
  \item{seed}{optional starting value for random generator. Default is seed = 1000.} 
}
\details{
  This function computes EMVE as described in Danilov et al. (2012). Two subsampling schemes can be used for computing EMVE:
  uniform subsampling and the clustering-based subsampling as described in Leung and Zamar (2016).
  For uniform subsampling, the number of subsamples must be large to ensure high breakdown point. For clustering-based 
  subsampling, the number of subsamples can be smaller. The subsample size \eqn{n_0} must be chosen to be larger than \eqn{p} to avoid singularity. 
  
  In the algorithm, there exists a concentration step in which Gaussian MLE is computed for \eqn{50\%} of the data points using the 
  classical EM-algorithm multiplied by a scalar factor. This step is repeated for each subsample. As the computation can be heavy as the 
  number of subsample increases, we set by default the maximum number of iteration of classical EM-algorithm (i.e. \code{maxits}) as 5. 
  Users are encouraged to refer to Danilov et al. (2012) for details about the algorithm
  and Rubin and Little (2002) for the classical EM-algorithm for missing data. 
}
\value{
  An S4 object of class \code{\link{emve-class}} which is a subclass of the virtual class \code{\link{CovRobMissSc-class}}.
  The output S4 object contains the following slots:
  
    \tabular{ll}{
    \code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
    \code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
    \code{sc} \tab Estimated EMVE scale. Can be accessed via \code{\link{getScale}}. \cr
    \code{pmd} \tab Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. \cr
    \code{pmd.adj} \tab Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. \cr
    \code{pu} \tab Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. \cr
    \code{call} \tab Object of class \code{"language"}. Not meant to be accessed. \cr
    \code{x} \tab Input data matrix. Not meant to be accessed. \cr
    \code{p} \tab Column dimension of input data matrix. Not meant to be accessed. \cr
    \code{estimator}\tab Character string of the name of the estimator used. Not meant to be accessed. \cr
    }
}
\references{
    Danilov, M., Yohai, V.J., Zamar, R.H. (2012). Robust Esimation of Multivariate
    Location and Scatter in the Presence of Missing Data. \emph{Journal of the American Statistical Association}
    \bold{107}, 1178--1186.
    
    Leung, A. and Zamar, R.H. (2016). Multivariate Location and Scatter Matrix Estimation Under Cellwise and Casewise Contamination.
    Submitted.
    
    Rubin, D.B. and Little, R.J.A. (2002). \emph{Statistical analysis with missing data} (2nd ed.). New York: Wiley.
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Ruben H. Zamar, Mike Danilov, Victor J. Yohai}
\seealso{
  \code{\link{GSE}}, \code{\link{emve-class}}
}
