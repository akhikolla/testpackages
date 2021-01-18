\name{law.cstr}
\alias{law.cstr}
\title{
Gives information about a given law.
}
\description{
To obtain the name of a law as well as its default number of parameters
and default parameter values.
}
\usage{
law.cstr(law.index, law.pars = NULL)
}
\arguments{
  \item{law.index }{
a single integer value corresponding to the index of a distribution as
given by function \code{\link{getindex}}.
}
  \item{law.pars }{
vector of the values of the parameters of the law specified in
\code{law.index}. If \code{NULL}, the default values are used.
}
}
\details{
This function can be useful to construct a title for a graph for example. 
}
\value{
\item{name }{name of the distribution with its parameters and the values
they take.}
\item{nbparams }{default number of parameters of the law.}
\item{law.pars }{values of the parameters.}
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{
See \code{\link{stat.cstr}}, \code{\link{getindex}},
\code{\link{getnbparlaws}}, \code{\link{getnbparstats}}.
}
\examples{
law.cstr(2)
}
\keyword{utilities}
