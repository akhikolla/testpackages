\name{RococoTestResults-class}
\docType{class}
\alias{RococoTestResults-class}
\alias{RococoTestResults}
\alias{rococotestresults}
\alias{show,RococoTestResults-method}

\title{Class "RococoTestResults"}
\description{S4 class for storing results of the robust rank correlation
test}
\section{Objects}{
Objects of this class can be created by calling \code{\link{rococo.test}}.
}
\section{Slots}{
The following slots are defined for \code{RococoTestResults} objects:
\describe{
    \item{\code{count}:}{number of times in which the test
      statistic for a random shuffle exceeded the test statistic of the
      true data; see \code{\link{rococo.test}}.}
    \item{\code{tnorm}:}{list identifying t-norm to use or
      two-argument function; see \code{\link{rococo}}.
      If one of the standard choices \code{"min"}, \code{"prod"}, or
      \code{"lukasiewicz"} has been used, the list has one component,
      \code{name} that contains the string identifying the t-norm.
      If a user-defined function has been used, the list has two
      components: \code{name} contains \code{"user-defined t-norm"}
      or the \code{name} attribute of the function object if available
      and \code{def} contains the function object itself.}
    \item{\code{input}:}{character string describing the input for which
      \code{\link{rococo.test}} has been called.}
    \item{\code{length}:}{number of samples for which
      \code{\link{rococo.test}} has been called.}
    \item{\code{p.value}:}{p-value of test.}
    \item{\code{p.value.approx}:}{p-value as based on a normal approximation of the null distribution.}
    \item{\code{r.values}:}{vector containing tolerance levels for the
      two inputs; see \code{\link{rococo.test}} or
      \code{\link{rococo}}.}
    \item{\code{numtests}:}{number of (random) shuffles performed by
      \code{\link{rococo.test}}.}
    \item{\code{exact}:}{logical indicating whether p-value has been
      computed exactly; see \code{\link{rococo.test}}.}
    \item{\code{similarity}:}{character (vector) identifying the
      similarity measure(s) used by \code{\link{rococo.test}}.}
    \item{\code{sample.gamma}:}{test statistic (robust gamma rank
      correlation coefficient) determined by
      \code{\link{rococo.test}}.}
    \item{\code{H0gamma.mu}:}{empirical mean of test
      statistic for random shuffles}
    \item{\code{H0gamma.sd}:}{empirical standard deviation of test
      statistic for random shuffles}
    \item{\code{perm.gamma}:}{in case \code{\link{rococo.test}} was
      called with \code{storeValues=TRUE}, this slot contains the
      vector of test statistics for random shuffles.}
    \item{\code{alternative}:}{alternative hypothesis used by
    \code{\link{rococo.test}}.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "RococoTestResults")}: d
      displays the most important information stored in
  \code{object}}
  }
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
\seealso{\code{\link{rococo.test}}, \code{\link{rococo}},
\code{\link{show-methods}}}
\examples{
## create data
f <- function(x) ifelse(x > 0.9, x - 0.9, ifelse(x < -0.9, x + 0.9, 0))
x <- rnorm(25)
y <- f(x) + rnorm(25, sd=0.1)

## perform correlation tests
ret <- rococo.test(x, y, similarity="classical", alternative="greater")
show(ret)

ret <- rococo.test(x, y, similarity="linear", alternative="greater")
show(ret)

ret <- rococo.test(x, y, similarity=c("classical", "gauss"),
                   r=c(0, 0.1), alternative="greater",
                   numtests=10000)
show(ret)
}
\keyword{htest}
\keyword{classes}
