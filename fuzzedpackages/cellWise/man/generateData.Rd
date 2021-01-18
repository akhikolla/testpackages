\name{generateData}
\alias{generateData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates artificial datasets with outliers
}
\description{
This function generates multivariate normal datasets with several possible types of outliers.
It is used in several simulation studies. For a detailed description, see the referenced papers.
}
\usage{
generateData(n, d, mu, Sigma, perout, gamma,
             outlierType = "casewise", seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{The number of observations
  }
  \item{d}{The dimension of the data.
  }
    \item{mu}{The center of the clean data.
  }
    \item{Sigma}{The covariance matrix of the clean data. Could be obtained from \code{\link{generateCorMat}}.
  }
  \item{outlierType}{The type of contamination to be generated. Should be one of:
  \itemize{
  \item \code{"casewise"}: Generates point contamination in the direction of the last eigenvector of \code{Sigma}.
  \item \code{"cellwisePlain"}: Generates cellwise contamination by randomly replacing a number of cells by \code{gamma}.
  \item \code{"cellwiseStructured"}: Generates cellwise contamination by first randomly sampling contaminated cells, after which for each row, they are replaced by a multiple of the smallest eigenvector of \code{Sigma} restricted to the dimensions of the contaminated cells.
  \item \code{"both"}: combines \code{"casewise"} and \code{"cellwiseStructured"}.
  }
  }
  \item{perout}{The percentage of generated outliers. For \code{outlierType = "casewise"} this is a fraction of rows. For \code{outlierType = "cellWisePlain"} or \code{outlierType = "cellWiseStructured"}, a fraction of \code{perout} cells are replaced by contaminated cells.
  For \code{outlierType = "both"}, a fraction of \eqn{0.5*}\code{perout} of rowwise
  outliers is generated, after which the remaining data is contaminated with a fraction of 
   \eqn{0.5*}\code{perout} outlying cells.
  }
  \item{gamma}{How far outliers are from the center of the distribution.
  }
 
  \item{seed}{
  Seed used to generate the data.
  }
  }
  
\value{
 A list with components:
  \itemize{
    \item{\code{X} \cr
    The generated data matrix of size \eqn{n \times d}.
}
\item{\code{indcells} \cr
A vector with the indices of the contaminated cells.
}
\item{\code{indrows} \cr
A vector with the indices of the rowwise outliers.
}
}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{generateCorMat}}
}

\examples{
n     <- 100
d     <- 5
mu    <- rep(0, d)
Sigma <- diag(d)
perout <- 0.1
gamma <- 10
data <- generateData(n, d, mu, Sigma, perout, gamma, outlierType = "cellwisePlain", seed  = 1)
pairs(data$X)
data$indcells
}