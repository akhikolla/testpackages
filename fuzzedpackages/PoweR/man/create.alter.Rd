\name{create.alter}
\alias{create.alter}
\title{
Create a list giving the type of test statistics.
}
\description{
Create a list giving the type of each test statistic for a given vector
of indices of these test statistics.
}
\usage{
create.alter(stat.indices = c(42, 51, 61), values.alter = NULL)
}
\arguments{
  \item{stat.indices }{vector of indices of test statistics, as given by
    function \code{\link{getindex}}.
}
  \item{values.alter }{
vector of the type of each test statistic in \code{stat.indices}. If
\code{NULL}, the default value will be returned.
}
}
\details{
See Section 3.3 in Lafaye de Micheaux, P. and Tran, V. A. (2014).
}
\value{
A named list. Each component of the list has the name of the
corresponding index in \code{stat.indices} (e.g. \code{stat}\eqn{xxx})
and has the value (in {0,1,2,3,4}) of the type of test (see Details above). 
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{
P. Lafaye de Micheaux, V. A. Tran
}
\seealso{
See \code{\link{getindex}}.
}
\examples{
create.alter()
}
\keyword{utilities}
