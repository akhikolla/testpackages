\name{compquant}
\alias{compquant}
\title{Computation of the quantile values only for one test statistic.}
\description{Functions for the computation of the quantile values only for one test statistic at a time and also one n value.}
\usage{
compquant(n,law.index,stat.index,probs=NULL,M=10^5,law.pars=NULL,
          stat.pars=NULL,model=NULL,Rlaw=NULL,Rstat=NULL,
center=FALSE, scale=FALSE)
}
\arguments{
  \item{n }{number of observations for each sample to be generated; length(\code{n})=1. This can also be set to 0 if you want to use
  your own function using the 'Rstat' argument (see below).}
  \item{law.index }{law index as given by \code{\link{getindex}}; length(\code{law.index})=1.}
  \item{stat.index }{stat index as given by \code{\link{getindex}}; length(\code{stat.index})=1.}
  \item{probs }{If not \code{NULL}, should be a vector of levels from which to compute the
    quantile values. If \code{NULL}, the levels 0.025,0.05,0.1,0.9,0.95,0.975
    will be used.}
  \item{M }{Number of Monte Carlo repetitions to use.}
  \item{law.pars }{\code{NULL} or a vector of length at most 4 containing 4 possible
  parameters to generate random values from distribution
  law(\eqn{\code{law.pars}[j],j<=4}). If \code{NULL}, the default parameter values for the law
  specified by \code{law.index} will be used.}
  \item{stat.pars }{A vector of parameters. 
    If NULL, the default parameter values for the statistic specified by this \code{stat.index} will be used.}
  \item{model }{NOT YET IMPLEMENTED. If \code{NULL}, no model is used.
  If an integer \eqn{i>0}, the model coded in the C function modele\eqn{i} is used.
  Else this should be an R function that
    takes three arguments: \code{eps} (vector of \eqn{\epsilon} values), \code{thetavec}
    (vector of \eqn{\theta} values) and \code{xvec} (vector or matrix of \eqn{x} values). This function should take a
    vector of errors, generate observations from a model (with
    parameters \code{thetavec} and values \code{xvec}) based on these errors, then
    compute and return the residuals from the model. See file 
    modele1.R in directory inst/doc/ for an example in multiple linear
    regression.}
  \item{Rlaw}{The user can provide its own (random generating) R
    function using this parameter. In this case, 'law.index' should be
    set to 0.}
  \item{Rstat }{If 'stat.index' is set to 0, an R function that outputs
    a list with components 'statistic' (value of the test statistic),
    'pvalue' (pvalue of the test; if not computable should be set to 0), 'decision' (1 if we reject the null,
    0 otherwise), 'alter' (see above), 'stat.pars' (see above),
    'pvalcomp' (1L if the pvalue can be computed, 0L otherwise),
    'nbparstat' (length of stat.pars).}
  \item{center}{Logical. Should we center the data generated}
  \item{scale}{Logical. Should we center the data generated}
}
\value{A list with \code{M} statistic values and also some quantiles (with
  levels 0.025,0.05,0.1,0.9,0.95,0.975), as well as the name of the law and the name of the test statistic used (just to be sure!).
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\examples{
compquant(n=50,law.index=2,stat.index=10,M=10^3)$quant
compquant(n=50,law.index=0,stat.index=10,M=10^3,Rlaw=rnorm)$quant}
\keyword{distribution}
\keyword{datagen}
