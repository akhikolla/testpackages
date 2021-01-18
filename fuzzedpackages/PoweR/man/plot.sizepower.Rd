\name{plot.sizepower}
\alias{plot.sizepower}
\title{
size-power curves.
}
\description{
This function produces a size-power curves plot.
}
\usage{
\method{plot}{sizepower}(x, xnull,legend.pos=NULL,...)
}
\arguments{
  \item{x }{
\code{Fx} object as returned by function \code{\link{calcFx}}.
}
  \item{xnull }{
\code{Fx} object as returned by function \code{\link{calcFx}}, but computed
under the null.
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
See \code{\link{plot.pvalue}}, \code{\link{plot.discrepancy}}, \code{\link{graph}}.
}
\examples{
## You can increase M for better results:
stind <- c(43,44,42)   # Indices of test statistics.
alter <-list(stat43=3,stat44=3,stat42=3) # Type for each test.
# Several p-values computed under the null.
pnull <- many.pval(stat.indices=stind,law.index=1,
                   n=100,M=100,N=10,alter=alter,null.dist=2,
                   method="MC")$pvals
Fxnull <- calcFx(pnull)
p <- many.pval(stat.indices=stind,law.index=4,n=100,
               M=100,N=10,alter=alter,null.dist=2,
               method="MC")$pvals
Fx <- calcFx(p)
plot.sizepower(Fx,Fxnull)
}
\keyword{hplot}
