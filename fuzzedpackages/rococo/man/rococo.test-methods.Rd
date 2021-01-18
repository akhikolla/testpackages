\name{rococo.test-methods}
\docType{methods}
\alias{rococo.test}
\alias{rococo.test-methods}
\alias{rococo.test,numeric,numeric-method}
\alias{rococo.test,formula,data.frame-method}
\title{Robust Gamma Rank Correlation Test}
\description{
  Methods performing a robust gamma rank correlation test
}
\usage{
\S4method{rococo.test}{numeric,numeric}(x, y,
            similarity=c("linear", "exp", "gauss",
                         "epstol", "classical"),
            tnorm="min", r=0, numtests=1000, storeValues=FALSE,
            exact=FALSE, alternative=c("two.sided", "less", "greater"),
            noVarReturnZero=TRUE)
\S4method{rococo.test}{formula,data.frame}(x, y, na.action, ...)
}
\arguments{
  \item{x}{a numeric vector or a formula; compulsory argument}
  \item{y}{compulsory argument; if \code{x} is a vector, \code{y}
    must be vector of the same length as \code{x}. If \code{x} is
    a formula, \code{y} must be a data frame.}
  \item{similarity}{a character string or a character vector identifying
    which type of similarity measure to use; see \code{\link{rococo}}
    for more details.}
  \item{tnorm}{t-norm used for aggregating results; see \code{\link{rococo}}
    for more details.}
  \item{r}{numeric vector defining the tolerances to be used;
    see \code{\link{rococo}} for more details.}
  \item{numtests}{number of random shuffles to perform; see details
    below.}
  \item{storeValues}{logical indicating whether the vector of test
    statistics should be stored in the output object (in slot
    \code{perm.gamma}, see \code{\linkS4class{RococoTestResults}}).}
  \item{exact}{logical indicating whether exact p-value should be
    computed; see details below.}
  \item{alternative}{indicates the alternative hypothesis and must be one of
    \code{"two.sided"}, \code{"greater"}, or
    \code{"less"}. Abbreviations are allowed as long as they are
    unique. \code{"greater"} corresponds to positive
    association, \code{"less"} to negative association.}
  \item{noVarReturnZero}{if \code{TRUE} (default), a correlation of 0
    and a p-value of 1 are returned if there is no variation in at least
    one of the two observables. Otherwise, a correlation of \code{NA}
    and a p-value of 1 are returned and a warning is issued.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}'s.  Defaults to \code{getOption("na.action")}.}
  \item{...}{all parameters specified are forwarded internally to the
    method \code{rococo.test} with signature \code{numeric,numeric}.}
}
\details{If called for numeric vectors, \code{rococo.test} computes the
  robust gamma rank correlation coefficient of \code{x} and \code{y}
  according to the specified parameters (see \code{\link{rococo}}) and
  then performs a permutation test to compute a p-value.
  If \code{exact=TRUE}, \code{rococo.test} attempts to compute an exact
  p-value and ignores the \code{numtests} argument. This is done by
  considering all possible permutations and computing the ratio of
  permutations for which the test statistic is at least as large/small
  as the test statistic for unshuffled data. This works only for 10 or
  less samples. Otherwise \code{exact=TRUE} is ignored, a warning is
  issued and random shuffles are considered to estimate the p-value (as
  follows next).
  If \code{exact=FALSE}, \code{numtests} random shuffles of \code{y} are
  performed and the empirical standard deviation of the robust gamma
  correlation values for these shuffled data sets is computed.
  Under the assumption that these values are normally distributed around
  mean zero, the p-value is then computed from this distribution in the
  usual way.
  Note that a too small choice of the number of shuffles (parameter
  \code{numtests}) leads to unreliable p-values.

  If \code{rococo.test} is called for a formula \code{x} and a data
  frame \code{y}, then the method checks whether the formula \code{x} correctly
  extracts two columns from \code{y} (see examples below). If so, the
  two columns are extracted and the robust gamma rank correlation test
  is applied to them according to the specified parameters.

  Note that \code{exact=TRUE} may result in long computation times for
  user-defined t-norms.
}
\value{
  Upon successful completion, the function returns an object of class
  \code{\linkS4class{RococoTestResults}} containing the results.
}
\author{Martin Krone & Ulrich Bodenhofer
\email{rococo@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/rococo/}

U. Bodenhofer, M. Krone, and F. Klawonn (2013).
Testing noisy numerical data for monotonic association.
\emph{Inform. Sci.} \bold{245}:21-37. DOI: 10.1016/j.ins.2012.11.026.

U. Bodenhofer and F. Klawonn (2008).
Robust rank correlation coefficients on the basis of fuzzy orderings:
initial steps. \emph{Mathware Soft Comput.} \bold{15}(1):5-20.}
\seealso{\code{\link{rococo}}}
\examples{
## create data
f <- function(x) ifelse(x > 0.9, x - 0.9, ifelse(x < -0.9, x + 0.9, 0))
x <- rnorm(25)
y <- f(x) + rnorm(25, sd=0.1)

## perform correlation tests
rococo.test(x, y, similarity="classical", alternative="greater")
rococo.test(x, y, similarity="linear", alternative="greater")
rococo.test(x, y, similarity=c("classical", "gauss"), r=c(0, 0.1),
            alternative="greater", numtests=10000)

## the formula variant
require(datasets)
data(iris)
rococo.test(~ Petal.Width + Petal.Length, iris, similarity="linear",
            alternative="two.sided")
}
\keyword{htest}
\keyword{methods}