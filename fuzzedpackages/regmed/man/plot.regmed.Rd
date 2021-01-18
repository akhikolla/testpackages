\name{plot.regmed}
\alias{plot.regmed}
\title{
Plot a regmed object
}
\description{
Plot regmed object with arrows from exposure (x) to mediators, and from
mediators to outcome (y). Mediators are plotted from top to bottom
according to the size of their effect (size =  |alpha * beta| )
}
\usage{
\method{plot}{regmed}(x, cex=.9, lwd=1.5, lty=1, \dots)
}
\arguments{
  \item{x}{ regmed object, returned from \code{regmed.fit()} }
  \item{cex}{ character size for labels in the diagram }
  \item{lwd}{ line width for arrows in the diagram }
  \item{lty}{ line type for arrows in the diagram }
  \item{\dots}{ optional arguments for plot method }
}
\details{
Arrows are colored red when corresponding alpha/beta is positive, and
blue when corresponding alpha/beta is negative. Same for delta, direct effect.
}
\value{nothing is returned}
\author{
Daniel Schaid
}
\seealso{
  \code{\link{regmed.fit}}
  \code{arrows}
}
\examples{
data(regmed_example)

y <- regmed_example$y
x <- regmed_example$x
med <- regmed_example[, -c(1,2)]
fit.grid <- regmed.grid(x, med, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)

fit.trim <- trim.best(fit.grid)
which.med <- colnames(med) %in% dimnames(fit.trim$alpha)[[1]]
med.selected <- med[, which.med]

fit.regmed <- regmed.fit(x, med.selected, y, lambda = 0.2, frac.lasso=.8)
plot(fit.regmed)
}