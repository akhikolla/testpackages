\name{powcomp.easy}
\alias{powcomp.easy}
\title{Computation of power and level tables for hypothesis tests.}
\description{Functions for the computation of power and level tables for hypothesis tests, in LaTeX format.}
\usage{
powcomp.easy(params,M=10^5,model=NULL,Rlaws=NULL,Rstats=NULL,center=FALSE, scale=FALSE)
}
\arguments{
  \item{M }{number of Monte Carlo repetitions to use.}
  \item{params }{
    matrix with (at least) 11 named-columns with names (\code{n}, \code{law}, \code{stat},
    \code{level}, \code{cL}, \code{cR}, \code{alter}, \code{par1}, \code{par2}, \code{par3}, \code{par4}). Each row of
    \code{params} gives the necessary parameters for a simulation of powers.

    \code{n} : sample size;
    
    \code{law} : integer giving the index of the law considered;
    
    \code{stat} : integer giving the index of the test statistic
    considered (can be 0 if you want to use your own function for some
    test statistic; see 'Rstats' argument);
    
    \code{level} : double, this is the significance level desired;
    
    \code{cL} : left critical value (can be NA);
    
    \code{cR} : right critical value (can be NA);
    
    \code{alter} : type of test (integer value in {0,1,2,3,4});
    
    \code{par}j: values of the parameters of the distribution specified by
    \code{law} (can be \code{NA}).
    
    See 'Details section'.
  }
\item{model }{NOT YET IMPLEMENTED. If \code{NULL}, no model is used.
  If an integer \eqn{i>0}, the model coded in the C function modele\eqn{i} is used.
  Else this should be an R function that
    takes three arguments: \code{eps} (vector of \eqn{\epsilon} values), \code{thetavec}
    (vector of \eqn{\theta} values) and \code{xvec} (vector or matrix of
  \eqn{x} values). This function should take a
    vector of errors, generate observations from a model (with
    parameters \code{thetavec} and values \code{xvec}) based on these errors, then
    compute and return the residuals from the model. See function
    modele1.R in directory inst/doc/ for an example in multiple linear
    regression}
  \item{Rlaws}{When some law indices in second column of 'params' are
  equal to 0, this means that you will be using some R random generators
  not hardcoded in C in the package. In that case, you should provide
  the names of the random generation functions in the corresponding
  components of a list; the other components should be set to NULL.}
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
\details{
  If both \code{cL} and \code{cR} are \code{NA}, no critical values are used
    and the decision to reject (or not) the hypothesis is taken using
    the \eqn{p}-value.

    If a test statistic depends upon some parameters, these can be added
    (in a correct order) in the last columns of \code{params}. If other
    test statistics are considered simultaneously (in the same
    \code{params} matrix) and if not all the test statistics have the
    same number of parameters, \code{NA} values should be used to
    complete empty cells of the matrix.
}
\value{The powers for the different statistics and laws specified in the rows
  of \code{params}, NOT YET provided in the form of a LaTeX table. This version is easier to use (but slower) 
  than the \code{\link{powcomp.fast}} version. It should be used in the process of investigating the power of test statistics under
  different alternatives. But when you are ready to produce results
  for publication in a paper, please use the \code{\link{powcomp.fast}} version and
  its \code{print} method..
}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\examples{
# Warning: the order of the parameters of the law (4 maximum) is important!
sim1 <- c(n=100,law=2,stat=10,level=0.05,cL=NA,cR=0.35,alter=3,
          par1= 2.0,par2=NA,par3=NA,par4=NA,parstat1=NA,parstat2=NA)
sim2 <- c(n=100,law=2,stat=17,level=0.10,cL=-0.30,cR=NA,alter=1,
          par1=-1.0,par2=3.0,par3=NA,par4=NA,parstat1=NA,parstat2=NA)
sim3 <- c(n=100,law=2,stat=31,level=0.10,cL=NA,cR=0.50,alter=3,
          par1=-1.0,par2=3.0,par3=NA,par4=NA,parstat1=0.7,parstat2=NA)
sim4 <- c(n=100,law=7,stat=80,level=0.10,cL=NA,cR=9.319,alter=3,
          par1=NA,par2=NA,par3=NA,par4=NA,parstat1=1,parstat2=5)
params <- rbind(sim1,sim2,sim3,sim4)
powcomp.easy(params,M=10^2)
sim5 <- c(n=100,law=0,stat=80,level=0.10,cL=NA,cR=9.319,alter=3,
          par1=NA,par2=NA,par3=NA,par4=NA,parstat1=1,parstat2=5)
params <- rbind(params,sim5)
powcomp.easy(params,M=10^2,Rlaws=list(NULL,NULL,NULL,NULL,rnorm))
}
\keyword{datagen}
\keyword{htest}
