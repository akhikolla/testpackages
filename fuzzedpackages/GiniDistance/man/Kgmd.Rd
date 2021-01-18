\name{Kgmd}
\title{ Kernel Gini Mean Difference Statistics}
\alias{Kgmd}
\description{
  Computes Kernel Gini mean difference statistics,
  in which Xs are quantitative, sigma is kernel standard deviation, alpha is an exponent on the Euclidean distance and returns the kernel Gini mean difference.
}
\usage{
  Kgmd(x, sigma)
}
\arguments{
  \item{x}{ data }
  \item{sigma}{ kernel standard deviation}
}
\details{
  \code{Kgmd} compute kernel Gini mean difference statistics for data.
  It is a self-contained R function dealing with both univariate and multivariate data. 
  
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Argument
  \code{x}, is treated as data.
  
  Energy distance based statistics naturally generalizes from a Euclidean space to metric spaces (\cite{Lyons13}). By using a positive definite kernel (Mercer kernel) (\cite{Mercer1909}), distributions are mapped into a RKHS (\cite{Smola07}) with a kernel induced distance. Hence one can extend energy distances to a much richer family of statistics defined in RKHS (\cite{Sejdinovic13}). Let \eqn{\kappa: R^q \times R^q \rightarrow R} be a Mercer kernel (\cite{Mercer1909}). There is an associated RKHS \eqn{H_{\kappa}} of real functions on \eqn{R^q} with reproducing kernel \eqn{\kappa}, where the function \eqn{d: {R}^q \times {R}^q \rightarrow{{R}}} defines a distance in \eqn{\mathcal{H}_\kappa},
  \deqn{
  d_\kappa(x,x') = \sqrt{\kappa(x,x) + \kappa(x',x') - 2 \kappa(x,x')}.}

Here \code{Kgcov} is defined as Gini distance covariance between \eqn{x} and \eqn{\mathrm{rank}(x)}. 

}
\value{
  \code{Kgmd} returns the sample Kernel Gini distance
}
\seealso{
  \code{\link{gCov}}  \code{\link{gCor}}  \code{\link{dCor}}
}
\references{
  
  Lyons, R. (2013). Distance covariance in metric spaces. The Annals of Probability, 41 (5), 3284-3305. 
  
  
}

\examples{
  x<-iris[,1]
  Kgmd(x, sigma=1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
