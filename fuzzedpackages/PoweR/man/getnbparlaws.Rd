\name{getnbparlaws}
\alias{getnbparlaws}
\title{
Retrieve the default number of parameters of some laws.
}
\description{
Retrieve the default number of parameters of the distributions in the package.
}
\usage{
getnbparlaws(law.indices = NULL)
}
\arguments{
  \item{law.indices }{vector of the indices of the distributions from which to retrieve the
default number of parameters. If \code{NULL}, all the distributions will
be considered.
}
}
\value{
  The default number of parameters for the laws specified in
  \code{law.indices}.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{
See \code{\link{getnbparstats}}, \code{\link{getindex}},
\code{\link{law.cstr}}, \code{\link{stat.cstr}}.
}
\examples{
## Default numbers of parameters for all the distributions in the package:
getnbparlaws()
## The Gaussian distribution has two parameters:
getnbparlaws(2)
}
\keyword{utilities}
