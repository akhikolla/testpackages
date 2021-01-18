\name{calcFx}
\alias{calcFx}
\title{
Empirical distribution function of p-values.
}
\description{
This function computes, at given points, the value of the empirical
distribution function of a sample of \eqn{p}-values.
}
\usage{
calcFx(pval.mat, x = c(seq(0.001, 0.009, by = 0.001), seq(0.01, 0.985,by = 0.005),
       seq(0.99, 0.999, by = 0.001)))
}
\arguments{
  \item{pval.mat }{matrix whose each column contains a vector of p-values for a given test
statistic. The column names of this matrix should be set to the names of
the various test statistics considered, whereas the rownames should all be
set to the name of the distribution under which the p-values have been
computed. This matrix can be obtained using function \code{\link{many.pval}}.
}
  \item{x }{vector of points at which to evaluate the empirical distribution function.
}
}
\details{
See equation (2) in Lafaye de Micheaux and Tran (2014).
}
\value{
An object of class \code{Fx} is returned, which contains a list whose
components are:
  \item{Fx.mat }{matrix whose ith column contains the values of the empirical distribution function (evaluated
  at the points in vector \code{x}) of the \eqn{p}-values of the ith test statistic.}
  \item{x }{same vector \code{x} as input.}
  \item{law }{name of the distribution under which the p-values have
    been computed. Should correspond to
    the row names of \code{pval.mat}.}
  \item{statnames }{names of the test statistics. Should correspond to
    the column names of \code{pval.mat}.}
  \item{N }{number of \eqn{p}-values used.}
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{
P. Lafaye de Micheaux, V. A. Tran
}
\seealso{
\code{\link{many.pval}}, \code{\link{plot.pvalue}}, \code{\link{plot.discrepancy}}, \code{\link{plot.sizepower}}
}
\examples{
stind <- c(43, 44, 42)   # Indices of test statistics.
alter <-list(stat43 = 3, stat44 = 3, stat42 = 3) # Type for each test.
# Several p-values computed under the null.
pnull <- many.pval(stat.indices = stind, law.index = 1,
                        n = 100, M = 10, N = 10, alter = alter,
                        null.dist = 1,
                        method = "direct")$pvals
xnull <- calcFx(pnull)
}
\keyword{nonparametric}
