\name{statcompute}
\alias{statcompute}
\title{Performs a hypothesis test for the given value of statistic.}
\description{Performs the hypothesis test for those added in the package.}
\usage{
statcompute(stat.index, data, levels = c(0.05,0.1), critvalL = NULL,
            critvalR = NULL, alter = 0, stat.pars = NULL, check = TRUE)
}
\arguments{
  \item{stat.index }{one statistic index as given by function
    \code{\link{getindex}}.}
  \item{data }{sample from which to compute the statistic.}
  \item{levels }{vector of desired significance levels for the test.}
  \item{critvalL }{\code{NULL} or vector of left critival values.}
  \item{critvalR }{\code{NULL} or vector of right critival values.}
  \item{alter }{0: two.sided=bilateral, 1: less=unilateral, 2:
  greater=unilateral, 3: bilateral test that rejects H0 only for large
  values of the test statistic, 4: bilateral test that rejects H0 only
  for small values of the test statistic.} 
  \item{stat.pars }{a vector of parameters. 
    If \code{NULL}, the default parameter values for this statistic will
    be used.}
  \item{check }{Logical. If \code{FALSE} it will execute much faster,
    but in this case be sure to give a value to the
    'stat.pars' argument; if you don't know what value to give, use \code{rep(0.0,
      getnbparstats(stat.index))} as a default value.}
 }
 \value{A list with components:
   \item{statistic }{the test statistic value}
   \item{pvalue }{the \eqn{p}-value}
   \item{decision }{the vector of decisions, same length as \code{levels}}
   \item{alter }{\code{alter}}
   \item{stat.pars }{\code{stat.pars}}
   \item{symbol }{how the test is noted}
 }
 \details{
The function statcompute() should not be used in simulations since it is
NOT fast. Consider instead using \code{powcomp.easy} or
\code{powcomp.fast}. Se also in the Example section below for a fast
approach using the \code{.C} function (but be warned that giving wrong
values of arguments can crash your session!).

   }
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\examples{
data <- rnorm(50)
# Shapiro-Wilk test:
statcompute(21, data, levels = c(0.05, 0.1), critvalL = NULL, critvalR = NULL,
            alter = 0, stat.pars = NULL)
# Identical to:
shapiro.test(data)

# The function statcompute() should not be used in simulations since it
#  is NOT fast. Consider instead the call below (but see the Details
#  Section):
.C("stat21", data = data, n = 50L, levels = 0.05, nblevels = 1L, name =
rep(" ", 50), getname = 0L, statistic = 0, pvalcomp = 1L, pvalue = 0, cL = 0.0,
cR = 0.0, usecrit = 0L, alter = 4L, decision = 0L, stat.pars = 0.0, 
nbparstat = 0L)

}
\keyword{htest}
