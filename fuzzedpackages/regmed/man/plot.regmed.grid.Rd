\name{plot.regmed.grid}
\alias{plot.regmed.grid}
\title{
Plots for regmed.grid object.
}
\description{
Creates 2 plots: (1) BIC vs. lambda, and (2) Coefficients Alpha and Beta of mediator vs. lambda.
}
\usage{
\method{plot}{regmed.grid}(x, as.log=FALSE, \dots)
}
\arguments{
  \item{x}{ regmed.grid object, returned by the \code{regmed.grid()} function }
  \item{as.log}{Logical; if TRUE,  plot lambda on the log scale}
  \item{\dots}{optional arguments for plot method}
}
\value{nothing is returned}
\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}
\seealso{
\code{\link{regmed.grid}},
}
\examples{
data(regmed_example)
y <- regmed_example$y
x <- regmed_example$x
med <- regmed_example[, -c(1,2)]
fit.grid <- regmed.grid(x, med, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)
plot(fit.grid)
}
