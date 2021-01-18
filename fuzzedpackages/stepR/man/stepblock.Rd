\name{stepblock}
\alias{stepblock}
\alias{[.stepblock}
\alias{print.stepblock}
\alias{plot.stepblock}
\alias{lines.stepblock}

\title{Step function}

\description{
Constructs an object containing a step function sampled over finitely many values.
}

\usage{
stepblock(value, leftEnd = c(1, rightEnd[-length(rightEnd)] + 1), rightEnd, x0 = 0)
\method{[}{stepblock}(x, i, j, drop = if(missing(i)) TRUE else if(missing(j)) FALSE else length(j) == 1, \dots)
\method{print}{stepblock}(x, \dots)
\method{plot}{stepblock}(x, type = "c", xlab = "x", ylab = "y", main = "Step function", sub = NULL, \dots)
\method{lines}{stepblock}(x, type = "c", \dots)
}

\arguments{
\item{value}{a numeric vector containing the fitted values for each block; its length gives the number of blocks}
\item{leftEnd}{a numeric vector of the same length as \code{value} containing the left end of each block}
\item{rightEnd}{a numeric vector of the same length as \code{value} containing the right end of each block}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started, i.e. before \code{leftEnd[1]}}
\item{x}{the object}
\item{i, j, drop}{see \code{\link{[.data.frame}}}
\item{type}{\code{"c"} to plot jumps in the middle between the end of the previous block (or \code{x0}) and the beginning of the following block; \code{"e"} to jump at the end of the previous block; \code{"b"} to jump at the beginning of the following block; capital letters also plot points}
\item{xlab, ylab, main, sub}{see \code{\link{plot.default}}}
\item{\dots}{for generic methods only}
}

\value{
For \code{stepblock} an object of class \code{stepblock}, i.e. a \code{\link{data.frame}} with columns \code{value}, \code{leftEnd} and \code{rightEnd} and \code{\link{attr}}ibute \code{x0}.
}

\note{
For the purposes of this package step functions are taken to be left-continuous, i.e. the function jumps \strong{after} the \code{rightEnd} of a block.

However, step functions are usually sampled at a discrete set of points so that the exact position of the jump is unknown, except that it has to occur before the next sampling point; this is expressed in the implementation by the specification of a \code{leftEnd} \strong{within} the block so that every \code{rightEnd} and \code{leftEnd} is a sampling point (or the boundary of the observation window), there is no sampling point between one block's \code{rightEnd} and the following block's \code{leftEnd}, while the step function is constant at least on the closed interval with boundary \code{leftEnd}, \code{rightEnd}.
}

\seealso{\code{\link{step}}, \code{\link{stepfit}}, \link{family}, \code{\link{[.data.frame}}, \code{\link{plot}}, \code{\link{lines}}}

\examples{
# step function consisting of 3 blocks: 1 on (0, 3]; 2 on (3, 6], 0 on (6, 8]
# sampled on the integers 1:10
f <- stepblock(value = c(1, 2, 0), rightEnd = c(3, 6, 8))
f
# show different plot types
plot(f, type = "C")
lines(f, type = "E", lty = 2, col = "red")
lines(f, type = "B", lty = 3, col = "blue")
legend("bottomleft", legend = c("C", "E", "B"), lty = 1:3, col = c("black", "red", "blue"))
}

\keyword{nonparametric}
