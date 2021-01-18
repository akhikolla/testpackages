\name{estLocScale}
\alias{estLocScale}

\title{
Estimate robust location and scale
}
\description{
Estimate a robust location estimate and scale estimate of every column in \code{X}.
}
\usage{
estLocScale(X, type = "wrap", precScale = 1e-12,
center = TRUE, alpha = 0.5, nLocScale = 25000, silent = FALSE)
}

\arguments{
  \item{X}{
  The input data. It must be an \eqn{n} by \eqn{d} matrix or a data frame.
}
  \item{type}{
 The type of estimators used. One of:
 \itemize{
  \item{\code{"1stepM"}: \cr The location is the 1-step M-estimator with the biweight psi function. The scale estimator is the 1-step M-estimator using a Huber rho function with \eqn{b = 2.5}.
 }
  \item{\code{"mcd"}: \cr the location is the weighted univariate MCD estimator with cutoff \cr
  \eqn{\sqrt(qchisq(0.975,1))}. The scale is the corresponding weighted univariate MCD estimator, with a correction factor to make it approximately unbiased at gaussian data.}
 \item{\code{"wrap"}: \cr Starting from the initial estimates corresponding to option \code{"mcd"}, the location is the 1-step M-estimator with the wrapping psi function with \eqn{b = 1.5} and \eqn{c = 4}. The scale estimator is the same as in option \code{"mcd"}.
 }

 }
 Defaults to "wrap".
}

  \item{precScale}{
The precision scale used throughout the algorithm. Defaults to \eqn{1e-12}.
}

\item{center}{
Whether or not the data has to be centered before calculating the scale. Not in use for \code{type = "mcd"}. Defaults to \code{TRUE}.
}

\item{alpha}{
The value of \eqn{\alpha} in the univariate mcd, must be between 0.5 and 1. The subsetsize is \eqn{h = \lceil \alpha n \rceil}. Only used for \code{type = "mcd"}. Defaults to \eqn{\alpha = 0.5}.
}
\item{nLocScale}{
 If \code{nLocScale} \eqn{< n}, \code{nLocScale} observations are sampled to compute the location and scale. This speeds up the computation if \eqn{n} is very large. When \code{nLocScale} \eqn{= 0} all observations are used. Defaults to \code{nLocScale} \eqn{= 25000}.
}
\item{silent}{
Whether or not a warning message should be printed when very small scales are found. Defauts to \code{FALSE}.
}
}


\value{
A list with components: \cr
  \itemize{
    \item{\code{loc} \cr
   A vector with the estimated locations.
    }
    \item{\code{scale} \cr
   A vector with the estimated scales.
    }
}
}

\references{
Raymaekers, J., Rousseeuw P.J. (2019). Fast robust correlation for high dimensional data. \emph{Technometrics}, published online. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2019.1677270}{(link to open access pdf)}
}
\author{
 Raymaekers, J. and Rousseeuw P.J.
}

\seealso{
\code{\link{wrap}}
}

\examples{
library(MASS) 
set.seed(12345) 
n = 100; d = 10
X = mvrnorm(n, rep(0, 10), diag(10))
locScale = estLocScale(X)
# For more examples, we refer to the vignette:
vignette("wrap_examples")
}