\name{pvalueMC}
\alias{pvalueMC}
\title{
Monte-Carlo computation of a p-value for one single test statistic.
}
\description{
This function can compute the p-value associated with a test statistic
value from a sample of observations.
}
\usage{
pvalueMC(data, stat.index, null.law.index, M = 10^5, alter, null.law.pars = NULL,
         stat.pars = NULL, list.stat = NULL, method = c("Fisher"),
         center = FALSE, scale = FALSE)
}
\arguments{
  \item{data }{sample of observations.}
  \item{stat.index }{index of a test statistic as given by function \code{\link{getindex}}.}
  \item{null.law.index }{index of the distribution to be tested (the
    null hypothesis distribution), 
    as given by function \code{\link{getindex}}.}
  \item{M }{number of Monte-Carlo repetitions to use.}
  \item{alter }{value (in {0,1,2,3,4}) giving the the type of test (See Section 3.3 in Lafaye de Micheaux, P. and Tran, V. A. (2014)). }
  \item{null.law.pars }{vector of parameters for the law. The length of this
    parameter should not exceed 4. If not provided, the default values
    are taken using \code{\link{getindex}} function.}
  \item{stat.pars }{a vector of parameters. 
    If \code{NULL}, the default parameter values for the statistic specified by this statistic wil be used.
}
  \item{list.stat }{if not \code{NULL}, a vector of test statistic values should be
provided. If \code{NULL}, these values will be computed.}
  \item{method }{method to use for the computation of the \eqn{p}-value. Only 'Fisher' method
is available for the moment.}
  \item{center}{Logical. Should we center the data generated}
  \item{scale}{Logical. Should we center the data generated}
}
\value{The Monte-Carlo p-value of the test.}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{statcompute}}.}
\examples{
x <- rnorm(100)
statcompute(1,x,level = c(0.05),alter = 3)$pvalue
pvalueMC(x,stat.index = 1,null.law.index = 2,M = 10^5,alter = 3)
}
\keyword{univar}
\keyword{htest}

