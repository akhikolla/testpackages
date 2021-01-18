\name{powcomp.fast}
\alias{powcomp.fast}
\title{Computation of power and level tables for hypothesis tests.}
\description{Functions for the computation of power and level tables for
  hypothesis tests, with possible use of a cluster.}
\usage{
powcomp.fast(law.indices,stat.indices,vectn = c(20,50,100),M = 10^3,levels = c(0.05,0.1),
         critval = NULL,alter = create.alter(stat.indices),parlaws = NULL, 
         parstats = NULL,nbclus = 1,model = NULL,null.law.index = 2,null.law.pars = NULL,
         Rlaws=NULL, Rstats = NULL, center=FALSE, scale=FALSE)
}
\arguments{
  \item{law.indices }{vector of law indices as given by function \code{\link{getindex}}.}
  \item{stat.indices }{vector of statistic indices as given by function \code{\link{getindex}} (some components can be 0 if you want to use your own function for some
    test statistics; see 'Rstats' argument).}
  \item{vectn }{vector of sample sizes (\eqn{n}) values.}
  \item{M }{number of Monte Carlo repetitions.}
  \item{levels }{vector of significance levels for the test.}
  \item{critval }{if not \code{NULL}, a named-list of critical values
    for each test statistic. The names of the list should be \code{stat}\eqn{j}, \eqn{j}
    taken in \code{stat.indices}. Note that if a single value of
    \code{critval$stat}\eqn{j} is povided, then it is the right critical value. If
    two values are provided, then these are the left and right critical
    values, in that order. If \code{NULL}, \code{critval} is computed
    using the function \code{\link{many.crit}}; in that case, be sure to
    provide the correct value for \code{null.law.index}.}
  \item{alter }{named-list of integer values (0: two.sided=bilateral, 1: less=unilateral, 2:
  greater=unilateral, 3: bilateral test that rejects H0 only for large
  values of the test statistic, 4: bilateral test that rejects H0 only
  for small values of the test statistic). The names of the list should be \code{stat}\eqn{j}, \eqn{j} taken in
  \code{stat.indices}.}
  \item{parlaws }{named-list of parameter values for each law to
    simulate. The names of the list should be \code{law}\eqn{j}, \eqn{j}
    taken in \code{law.indices}. The length of vector \code{law}\eqn{j} should not be greater than
    4 (we supposed than no common distribution has more than 4 parameters!).}
  \item{parstats }{named-list of parameter values for each statistic to simulate. 
    The names of the list should be \code{stat}\eqn{j}, \eqn{j} taken in \code{stat.indices} (in
    the same order). If \code{NULL}, the default parameter values for these statistics will be used.}
  \item{nbclus }{number of slaves to use for the computation on a
    cluster. This needs parallel or Rmpi package to be installed
    and functionnal on the system. Also the mpd daemon sould be started.}
  \item{model }{NOT YET IMPLEMENTED. If \code{NULL}, no model is used.
    If an integer \eqn{i>0}, the model coded in the C function modele\eqn{i} is used.
    Else this shoud be an R function that
    takes three arguments: \code{eps} (vector of \eqn{\epsilon} values), \code{thetavec}
    (vector of \eqn{\theta} values) and \code{xvec} (vector or matrix of
    \eqn{x} values). This function should take a
    vector of errors, generate observations from a model (with
    parameters \code{thetavec} and values \code{xvec}) based on these errors, then
    compute and return the residuals from the model. See function
    modele1.R in directory inst/doc/ for an example in multiple linear
    regression.}
  \item{null.law.index }{index of the law under the null. Only used, by
    \code{\link{many.crit}} function, if \code{critval} is \code{NULL}.}
  \item{null.law.pars }{vector of parameters corresponding to
    \code{null.law.index}.}
  \item{Rlaws }{When some law indices in 'law.indices' are equal to 0, this means that you will be using some R random generators. In that case, you should provide the names of the random generation functions in the corresponding components of 'Rlaws' list, the other components should be set to NULL.}
\item{Rstats}{A list. If in a given row of the 'params' matrix, the
  value of 'stat' is set to 0, the corresponding component of the list
  'Rstats' should be an R function that outputs
    a list with components 'statistic' (value of the test statistic),
    'pvalue' (pvalue of the test; if not computable should be set to 0), 'decision' (1 if we reject the null,
    0 otherwise), 'alter' (see above), 'stat.pars' (see above),
    'pvalcomp' (1L if the pvalue can be computed, 0L otherwise),
    'nbparstat' (length of stat.pars). If the value of 'stat' is not 0,
  then the corresponding component of 'Rstats' should be set to 'NULL'.}
  \item{center}{Logical. Should we center the data generated}
  \item{scale}{Logical. Should we center the data generated}
}
\value{A list of class \code{power} whose components are described below:
  \item{M }{number of Monte Carlo repetitions.}
  \item{law.indices }{vector of law indices as given by function \code{\link{getindex}}.}
  \item{vectn }{vector of sample sizes.}
  \item{stat.indices }{vector of test statistic indices as given by function \code{\link{getindex}}.}
  \item{decision }{a vector of counts (between 0 and \code{M}) of the
    decisions taken for each one of the \code{levels.len} * \code{laws.len} *
    \code{vectn.len} * \code{stats.len} combinations of (level,law,sample size,test statistic), to be
    understood in the following sense. The decision for the \eqn{l}-th level
    (in \code{levels}), \eqn{d}-th law (in \code{law.indices}), \eqn{n}-th sample
    size (in \code{vectn}) and \eqn{s}-th test statistic (in
    \code{stat.indices}) is given by:
    
    \code{decision}[\eqn{s} + \code{stats.len}*(\eqn{n}-1) +
                     \code{stats.len}*\code{vectn.len}*(\eqn{d}-1) + 
                      \code{stats.len}*\code{vectn.len}*\code{laws.len}*(\eqn{l}-1)]

    where \code{stats.len}, \code{vectn.len}, \code{laws.len} and \code{levels.len} are respectively the lengths of the vectors
    \code{stat.indices}, \code{vectn}, \code{law.indices} and \code{levels}.
  }
  \item{levels }{vector of levels for the test.}
  \item{cL }{left critical values used.}
  \item{cR }{right critical values used.}
  \item{usecrit }{a vector of 1s and 0s depending if a critical value
    has been used or not.}
  \item{alter }{type of each one of the tests in \code{stat.indices} used (0: two.sided=bilateral, 1: less=unilateral, 2:
          greater=unilateral, 3: bilateral test that rejects H0 only
          for large values of the test statistic, 4: bilateral test
          that rejects H0 only for small values of the test statistic).}
  \item{nbparlaws }{default number of parameters used for each law in \code{law.indices}.}
  \item{parlaws }{default values of the parameters for each law.}
  \item{nbparstats }{default number of parameters for each test
    statistic in \code{stat.indices}.}
  \item{parstats }{default values of the parameters for each test statistic.}
  \item{nbclus }{number of CPUs used for the simulations.}
}
\details{
This version is faster (but maybe less easy to
  use in the process of investigating the power of test statistics under
  different alternatives) than the \code{\link{powcomp.easy}} version.
 }
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\examples{
## Regenerate Table 6 from Puig (2000) (page 424)

law.index <- 1
# Take M = 50000 for accurate results 
M <- 10
vectn <- c(10,15,20,35,50,75,100)
level <- c(0.05)
stat.indices <- c(43,44,42,45,46)
law.indices <- c(2,3,4)
alter <- list(stat43 = 3,stat44 = 3,stat42 = 3,stat45 = 3,stat46 = 3)
critval <- many.crit(law.index,stat.indices,M,vectn,level,alter,
                     law.pars = NULL,parstats = NULL)
table6 <- powcomp.fast(law.indices,stat.indices,vectn,M,level,critval = critval,alter,
                       parlaws = NULL,parstats = NULL,nbclus = 1)
table6
}
\keyword{datagen}
\keyword{htest}
