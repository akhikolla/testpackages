\name{rococo}
\alias{rococo}
\title{Robust Gamma Rank Correlation Coefficient}
\description{
   Compute the robust gamma rank correlation coefficient
}
\usage{
rococo(x, y,
       similarity=c("linear", "exp", "gauss", "epstol", "classical"),
       tnorm="min", r=0, noVarReturnZero=TRUE)
}
\arguments{
  \item{x}{a numeric vector; compulsory argument}
  \item{y}{a numeric vector; compulsory argument; \code{x} and
    \code{y} need to have the same length}
  \item{similarity}{a character string or a character vector identifying
    which type of similarity measure to use; valid values are
    \code{"linear"} (default), \code{"exp"},
    \code{"gauss"}, \code{"epstol"}, and \code{"classical"}
    (abbreviations are allowed as
    long as they are unique). If \code{similarity} is
    a single string, the same similarity measure is taken for \code{x}
    and \code{y}. Different similarity measures can be used for \code{x}
    and \code{y} by supplying different similarity measures in
    \code{similarity[1]} and \code{similarity[2]}. Longer character
    vectors are allowed, but all but the first two entries are ignored.}
  \item{tnorm}{can be any of the following strings identifying a
    standard tnorm: \code{"min"} (minimum t-norm; default), \code{"prod"}
    (product t-norm), or \code{lukasiewicz} (Lukasiewicz t-norm);
    abbreviations are allowed as long as they are unique. Alternatively,
    \code{tnorm} can be a two-argument function defining a t-norm.}
  \item{r}{numeric vector defining the tolerances to be used; if a
    single value is supplied, the same value is used both for \code{x} and
    \code{y}. If a vector is supplied, \code{r[1]} is used as tolerance
    for \code{x} and \code{r[2]} is used as tolerance
    for \code{y}. If the classical crisp similarity is
    used, the corresponding entry/entries in \code{r} is/are ignored.
    Negative values are not allowed. Zeroes have a special
    meaning: if an entry in \code{r} is \code{0}, then the tolerance is
    automatically adapted to 10 percent of the interquartile range of
    the data.}
  \item{noVarReturnZero}{if \code{TRUE} (default), a correlation of 0
    is returned if there is no variation in at least one of the two
    observables. Otherwise, \code{NA} is returned and a warning is issued.}
}
\details{\code{rococo} computes the robust gamma rank correlation
  coefficient of \code{x} and \code{y} according to the specified
  parameters (see literature for more details).

  Note that \code{rococo} only works for \code{x} and \code{y} being
  numeric vectors, unlike the classical correlation measures implemented
  in \code{\link{cor}} which can also be computed for matrices or data
  frames.
}
\value{
  Upon successful completion, the function returns the robust gamma rank
  correlation coefficient.
}
\author{Martin Krone & Ulrich Bodenhofer
\email{rococo@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/rococo/}
  
U. Bodenhofer and F. Klawonn (2008).
Robust rank correlation coefficients on the basis of fuzzy orderings:
initial steps. \emph{Mathware Soft Comput.} \bold{15}(1):5-20.

U. Bodenhofer, M. Krone, and F. Klawonn (2013).
Testing noisy numerical data for monotonic association.
\emph{Inform. Sci.} \bold{245}:21-37. DOI: 10.1016/j.ins.2012.11.026.}
\seealso{\code{\link{rococo.test}}}
\examples{
## create data
f <- function(x) ifelse(x > 0.9, x - 0.9, ifelse(x < -0.9, x + 0.9, 0))
x <- rnorm(25)
y <- f(x) + rnorm(25, sd=0.1)

## compute correlation
rococo(x, y, similarity="classical")
rococo(x, y, similarity="linear")
rococo(x, y, similarity=c("classical", "gauss"), r=c(0, 0.1))
}
\keyword{htest}
