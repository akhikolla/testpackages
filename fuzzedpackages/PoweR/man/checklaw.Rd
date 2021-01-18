\name{checklaw}
\alias{checklaw}
\title{Check proper behaviour of a random generator}
\description{It is desirable to check if a newly added random generator coded in C
  behaves correctly. To perform this operation, one can superimpose the theoretical density on a histogram of the
  generated values.}
\usage{
checklaw(law.index, sample.size = 50000, law.pars = NULL, density =
NULL, trunc = c(-Inf, Inf), center = FALSE, scale = FALSE)
}
\arguments{
  \item{law.index }{index of the desired law, as given by \code{\link{getindex}}.}
  \item{sample.size }{number of observations to generate.}
  \item{law.pars }{vector of parameters for the law. The length of this
    parameter should not exceed 4. If not provided, the default values
    are used by means of \code{\link{getindex}} function.}
  \item{density }{a function of two arguments \code{x} and
    \code{pars}. Can also be a function of
    the arguments \code{x} and \code{pars[1]}, ..., \code{pars[k]}. See the two examples below.}
  \item{trunc }{vector of left and right truncation thresholds for the
    generated sample values. Only those values in between will be kept
    to build the histogram. This can be useful for a distribution
    with extreme values.}
  \item{center}{Logical. Should we center the data.}
\item{scale}{Should we scale the data.}
}
\value{Returns invisibly the data generated and make a plot showing histogram and density superimposed.
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran
}
\seealso{\code{\link{gensample}}}
\examples{
dlaplace1 <- function(x, mu, b) {dexp(abs(x - mu), 1 / b) / 2}
checklaw(1, density = dlaplace1)
dlaplace2 <- function(x, pars) {dexp(abs(x - pars[1]), 1 / pars[2]) / 2}
checklaw(1, density = dlaplace2)


checklaw(law.index = 2, sample.size = 50000, law.pars = c(2, 3), density
= dnorm)

## We use the 'trunc' argument to display the density in a region where
## no extreme values are present.
checklaw(27, density = dlaw27, trunc = c(-Inf,10))

# This one (Tukey) does not have a closed form expression for
# the density. But we can use the stats::density() function as
# follows.
res <- checklaw(18)
lines(density(res$sample), col = "blue")
}
\keyword{utilities}
\keyword{distribution}
