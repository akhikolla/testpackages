\name{wrap}
\alias{wrap}
\title{
Wrap the data.
}
\description{
Transforms multivariate data \code{X} using the wrapping function with \code{b = 1.5} and \code{c = 4}. By default, it starts by calling \code{\link{checkDataSet}} to clean the data and \code{\link{estLocScale}} to estimate the location and scale of the variables in the cleaned data. Alternatively, it works with user-provided vectors of location and scale given by \code{locX} and \code{scaleX}. 
}
\usage{
wrap(X, locX = NULL, scaleX = NULL, precScale = 1e-12,
     imputeNA = TRUE, checkPars = list())
}

\arguments{
  \item{X}{the input data. It must be an \eqn{n} by \eqn{d} matrix or a data frame.
}
  \item{locX}{The location estimates of the columns of the input data \code{X}. Must be a vector of length \eqn{d}.
}
  \item{scaleX}{The scale estimates of the columns of the input data \code{X}. Must be a vector of length \eqn{d}.
}
  \item{precScale}{The precision scale used throughout the algorithm. Defaults to \eqn{1e-12}
}
\item{imputeNA}{Whether or not to impute the \code{NA}s with the location estimate
of the corresponding variable. Defaults to \code{TRUE}.
}
\item{checkPars}{Optional list of parameters used in the call to
 \code{\link{checkDataSet}}. The options are:
  \itemize{
  \item{ \code{coreOnly} \cr
      If \code{TRUE}, skip the execution of checkDataset. Defaults to \code{FALSE}}
       \item{\code{numDiscrete}\cr
  A column that takes on numDiscrete or fewer values
  will be considered discrete and not retained in the cleaned data.
  Defaults to \eqn{5}.
}
  \item{\code{precScale} \cr
  Only consider columns whose scale is larger than precScale.
  Here scale is measured by the median absolute deviation.
  Defaults to \eqn{1e-12}.
}
 \item{\code{silent}\cr
  Whether or not the function progress messages should be printed.
  Defaults to \code{FALSE}.
}
      }
   }
}

\value{A list with components: \cr
  \itemize{
    \item{\code{Xw} \cr
   The wrapped data.
    }
    \item{\code{colInWrap} \cr 
    The column numbers of the variables which were wrapped. Variables which were filtered out by \code{\link{checkDataSet}} (because of a (near) zero scale for example), will not appear in this output.
    }
     \item{\code{loc} \cr 
    The location estimates for all variables used for wrapping.
    }
     \item{\code{scale} \cr 
    The scale estimates for all variables used for wrapping.
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
\code{\link{estLocScale}}
}
\examples{
library(MASS) 
set.seed(12345) 
n <- 100; d <- 10
X <- mvrnorm(n, rep(0, 10), diag(10))
locScale <- estLocScale(X)
Xw <- wrap(X, locScale$loc, locScale$scale)$Xw
# For more examples, we refer to the vignette:
vignette("wrap_examples")
}
