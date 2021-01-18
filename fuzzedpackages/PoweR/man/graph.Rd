\name{graph}
\alias{graph}
\title{
p-value plot, p-value discrepancy plot and size-power curves.
}
\description{
This function draws a \eqn{p}-value plot, a \eqn{p}-value discrepancy plot or a
size-power curves plot.
}
\usage{
graph(matrix.pval, xi = c(seq(0.001, 0.009, by = 0.001),
      seq(0.01, 0.985, by = 0.005), seq(0.99, 0.999, by = 0.001)),
      type = c("pvalue.plot", "pvalue.discrepancy", "size.power"),
      center = FALSE, scale = FALSE)
}
\arguments{
  \item{matrix.pval }{
a matrix of \eqn{p}-values as returned by function \code{\link{many.pval}}.
}
  \item{xi }{
a vector of values at which to compute the empirical distribution of
the \eqn{p}-values.
}
  \item{type }{
character. Indicate the type of plot desired.
}
  \item{center}{Logical. Should we center the data generated}
  \item{scale}{Logical. Should we center the data generated}
}
\details{
See Section 2.3 in Lafaye de Micheaux, P. and Tran, V. A. (2014).
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{
See \code{\link{plot.pvalue}}, \code{\link{plot.discrepancy}}, \code{\link{plot.sizepower}}.
}
\examples{
stind <- c(43,44,42)   # Indices of test statistics.
alter <-list(stat43=3,stat44=3,stat42=3) # Type for each test.
# Several p-values computed under the null.
# You can increase the values of M and N for better results.
matrix.pval <- many.pval(stat.indices=stind,law.index=1,
                        n=100,M=10,N=10,alter=alter,null.dist=1,
                        method="direct")
graph(matrix.pval)
}
\keyword{dplot}
\keyword{hplot}

