\name{many.pval}
\alias{many.pval}
\title{
Computes several \eqn{p}-values for many test statistics.
}
\description{
This function generates a sample of \code{n} observations from the law
specified in \code{law.index}. It then computes the value of each test
statistic specified in \code{stat.indices} and use it to obtain the
corresponding \eqn{p}-value under the null. The computation of these \eqn{p}-values
can be done using a Monte-Carlo simulation.
}
\usage{
many.pval(stat.indices, law.index, n = 100, M = 10^5, N = 100,
          alter = create.alter(stat.indices), law.pars = NULL, parstats = NULL,
          null.dist = 2, null.pars = NULL, method = c("direct", "MC"), Rlaw.index = NULL,
Rnull.dist = NULL, Rstats = NULL, center=FALSE, scale=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat.indices }{vector of test statistic indices as given by
    function \code{\link{getindex}} (some components can be 0 if you want to use your own function for some
    test statistics; see 'Rstats' argument).
}
  \item{law.index }{index of the distribution from which to generate
    observations used to compute the values of the test statistics specified with \code{stat.indices}.
}
  \item{n }{
integer. Size of the samples from which to compute the value of the test statistics.
}
  \item{M }{
integer. Number of Monte-Carlo repetitions. Only used when \code{method} = 'MC'.
}
  \item{N }{
integer. Number of \eqn{p}-values to compute for each test statistic.
}
  \item{alter }{
integer value in {0,1,2,3,4}. Type of test. See function \code{\link{create.alter}}.
}
  \item{law.pars }{
vector of the parameter values for the law specified in \code{law.index}.
}
  \item{parstats }{
named-list of vectors of parameters for the test statistics specified in
\code{stat.indices}. The names of the list should be \code{stat}\eqn{xxx} where \eqn{xxx}
are the indices specified in \code{stat.indices}.
}
  \item{null.dist }{
    used only if \code{method} = 'MC'. Integer value (as given by
    function \code{\link{getindex}}) specifying the distribution under
    the null hypothesis.
  }
    \item{null.pars }{vector of parameters for the null distribution}
  \item{method }{
    character. Either 'direct' to compute the \eqn{p}-value under the null using
for example the asymptotic distribution of the test statistic under the
null. This is not possible for all test statistics; or 'MC' to use a
Monte-Carlo simulation to approximate the distribution of the test
statistic under the null (specified by \code{null.dist}).
}
\item{Rlaw.index }{If 'law.index' is set to 0 then 'Rlaw.index' should be a (random generating) function.}
\item{Rnull.dist }{If 'null.dist' is set to 0 then 'Rnull.dist' should be a (random generating) function.}
\item{Rstats}{A list of same length as \code{stat.indices}. If a value of the vector \code{stat.indices} is set to 0, the corresponding component of the list
  \code{Rstats} should be an R function that outputs
    a list with components \code{statistic} (value of the test statistic),
    \code{pvalue} (pvalue of the test; if not computable should be set to 0), \code{decision} (1 if we reject the null,
    0 otherwise), \code{alter} (see above), \code{stat.pars} (see above),
    \code{pvalcomp} (1L if the pvalue can be computed, 0L otherwise),
    \code{nbparstat} (length of stat.pars). If a value of \code{stat.indices} is not 0,
  then the corresponding component of \code{Rstats} should be set to \code{NULL}.}
  \item{center}{Logical. Should we center the data generated}
  \item{scale}{Logical. Should we center the data generated}
}
\value{
\item{pvals }{the \code{N} x length(\code{stat.indices}) matrix of \eqn{p}-values.}
\item{stat.indices }{same as input.}
\item{n }{same as input.}
\item{M }{same as input.}
\item{alter }{same as input.}
\item{parstats }{same as input.}
\item{null.dist }{same as input.}
\item{method }{same as input.}
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{
See \code{\link{calcFx}}, \code{\link{graph}}.
}
\examples{
stind <- c(43,44,42)   # Indices of test statistics.
alter <-list(stat43=3,stat44=3,stat42=3) # Type for each test.
# Several p-values computed under the null.
# You can increase the values of M and N for better results.
matrix.pval <- many.pval(stat.indices=stind,law.index=1,
                        n=100,M=10,N=10,alter=alter,null.dist=1,
                        method="direct")

}
\keyword{htest}
