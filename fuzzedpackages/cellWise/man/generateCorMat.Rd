\name{generateCorMat}
\alias{generateCorMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates correlation matrices
}
\description{
This function generates correlation matrices frequently used in simulation
studies.
}
\usage{
generateCorMat(d, corrType = "ALYZ", CN = 100, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{The dimension of the correlation matrix. The resulting matrix is \eqn{d \times d}.
  }
  \item{corrType}{The type of correlation matrix to be generated. Should be one of:
  \itemize{
  \item \code{"ALYZ"}: Generates a correlation matrix as in Agostinelli et. al (2015).
  \item \code{"A09"}: Generates the correlation matrix defined by \eqn{\rho_{jh} = (-0.9)^{|h-j|}}.
  }
  Note that the option \code{"ALYZ"} produces a randomly generated correlation matrix.
  }
  \item{CN}{Condition number of the correlation matrix. Only used for \code{corrType = "ALYZ"}.
}
  \item{seed}{
  Seed used in \code{set.seed} before generating the correlation matrix. Only relevant for \code{corrType = "ALYZ"}.
  }
  }
  
\value{
A \eqn{d \times d} correlation matrix of the given type.
}
\references{
C. Agostinelli, Leung, A., Yohai, V. J., and Zamar, R. H. (2015).
Robust Estimation of Multivariate Location and Scatter in the Presence of Cellwise and Casewise Contamination. \emph{Test}, 24, 441-461.

Rousseeuw, P.J., Van den Bossche W. (2018). Detecting Deviating Data Cells. \emph{Technometrics}, \bold{60}(2), 135-145. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2017.1340909}{(link to open access pdf)}

J. Raymaekers and P.J. Rousseeuw (2020). Handling cellwise outliers by sparse
regression and robust covariance. \emph{Arxiv: 1912.12446}. \href{https://arxiv.org/abs/1912.12446}{(link to open access pdf)}
}
\author{
J. Raymaekers and P.J. Rousseeuw
}

\seealso{
\code{\link{generateData}}
}

\examples{
d     <- 5
Sigma <- generateCorMat(d, corrType = "ALYZ", seed = 1)
Sigma
}