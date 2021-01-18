\name{gauss.cor}
\alias{gauss.cor}
\title{Gauss Rank Correlation Estimator}
\description{
   Compute the Gaussian rank correlation estimate
}
\usage{gauss.cor(x, y)
}
\arguments{
  \item{x}{a numeric vector; compulsory argument}
  \item{y}{a numeric vector; compulsory argument; \code{x} and
    \code{y} need to have the same length}
}
\details{\code{gauss.cor} computes the Gaussian rank correlation estimate for
  \code{x} and \code{y}.

  Note that \code{gauss.cor} only works for \code{x} and \code{y} being
  numeric vectors, unlike the classical correlation measures implemented
  in \code{\link{cor}} which can also be computed for matrices or data
  frames.
}
\value{
  Upon successful completion, the function returns the Gaussian rank correlation estimate.
}
\author{Ulrich Bodenhofer
\email{rococo@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/rococo/}

K. Boudt, J. Cornelissen and C. Croux (2012).
The Gaussian rank correlation estimator: robustness properties.
\emph{Stat. Comput.} \bold{22}(2):471-483.}
\seealso{\code{\link{gauss.cor.test}}}
\examples{
## create data
f <- function(x) ifelse(x > 0.9, x - 0.9, ifelse(x < -0.9, x + 0.9, 0))
x <- rnorm(25)
y <- f(x) + rnorm(25, sd=0.1)

## compute correlation
gauss.cor(x, y)
}
\keyword{htest}
