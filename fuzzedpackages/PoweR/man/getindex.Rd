\name{getindex}
\alias{getindex}
\alias{print.index}
\title{Get indices of laws and statistics functions.}
\description{Print two correspondence tables between indices and random
  generators functions or test statistics functions programmed in C in this package.
  The first table gives indices/laws and the second one gives
  indices/statistics. These indices can be used in the 
  functions \code{\link{powcomp.easy}}, \code{\link{powcomp.fast}}, \code{\link{compquant}}, \code{\link{gensample}}, \code{\link{statcompute}}, \code{\link{checklaw}}.}
\usage{
getindex(law.indices = NULL, stat.indices = NULL)
}
\arguments{
  \item{law.indices }{if not \code{NULL}, select only the laws corresponding to
    this vector of indices.}
  \item{stat.indices }{if not \code{NULL}, select only the stats corresponding to
    this vector of indices.}
}
\value{A list with two matrices. The first one gives the correspondence
  between the indices and the laws (with also the number of parameters for
  each law as well as the default values). The second one gives the correspondence
  between the indices and the test statistics. Note that you can use the
  \code{law.indices} or \code{stat.indices} parameters of this function to obtain only
  some part of these tables of correspondence.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{getnbparlaws}}, \code{\link{getnbparstats}}, \code{\link{stat.cstr}},\code{\link{law.cstr}}.}
\examples{
getindex(1,c(4,3))
}
\keyword{utilities}
