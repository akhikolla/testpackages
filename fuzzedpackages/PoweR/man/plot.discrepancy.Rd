\name{plot.discrepancy}
\alias{plot.discrepancy}
\title{
\eqn{p}-value discrepancy plot.
}
\description{
This function produces a \eqn{p}-value discrepancy plot.
}
\usage{
\method{plot}{discrepancy}(x,legend.pos=NULL,...)
}
\arguments{
  \item{x }{
\code{Fx} object as returned by function \code{\link{calcFx}}.
}
\item{legend.pos }{
If \code{NULL}, position of the legend will be computed
automatically. Otherwise, it should be either a character vector in
"bottomright", "bottom", "bottomleft", "left", "topleft", "top",
"topright", "right" and "center". Or a numeric vector of length 2 giving
the x-y coordinates of the legend.
  }
\item{... }{further arguments passed to the \code{plot} or \code{points} functions.
  }
}
\details{
See Section 2.3 in Lafaye de Micheaux, P. and Tran, V. A. (2014).
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{
See \code{\link{plot.pvalue}}, \code{\link{plot.sizepower}}, \code{\link{graph}}.
}
\examples{
stind <- c(43,44,42)   # Indices of test statistics.
alter <-list(stat43=3,stat44=3,stat42=3) # Type for each test.
# Several p-values computed under the null.
pnull <- many.pval(stat.indices=stind,law.index=1,
                   n=100,N=10,alter=alter,null.dist=1,
                   method="direct")$pvals
xnull <- calcFx(pnull)
plot.discrepancy(xnull)
}
\keyword{hplot}
