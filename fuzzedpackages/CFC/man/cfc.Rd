\name{cfc}
\alias{cfc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cause-specific competing-risk survival analysis
}
\description{
Using adaptive generalized Newton-Cotes for calculating cumulative incidence functions.
}
\usage{
cfc(f.list, args.list, n, tout, Nmax = 100L, rel.tol = 1e-05, ncores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f.list}{In \code{R} mode, this is a list of survival functions, one per cause. Each survival function must have the prototype \code{f(t, args, n)}, where \code{t} is a vector of time-from-index values, \code{args} is a list of all arguments needed by the function, and \code{n} is the iterator that allows the function to produce a different output for each observation and/or Bayesian sample. The output is a vector of survival probabilities at times \code{t} for that particular observation/sample combination with iterator value \code{n}. In \code{C++} mode, this is a list of lists, one per cause. Each list must contain pointers to three \code{C++} functions, in order: 1) survival function of type \code{func} with prototype defined as \code{typedef arma::vec (*func)(arma::vec x, void* arg, int n)} (using \code{RcppArmadillo}'s \code{vec} class) and a similar interpetation as its \code{R} counterpart, 2) initializer function of type \code{initfunc}, with prototype defined as \code{typedef void* (*initfunc)(List arg)}, where \code{List} is the \code{Rcpp} wrapper class for \code{R} lists, and 3) resource de-allocator function of type \code{freefunc} with this prototype: \code{typedef void (*freefunc)(void *arg)}. See vignette for \code{C++} example.}
  \item{args.list}{List of arguments (each one a list), one per cause, to be supplied to the survival functions in \code{f.list}.}
  \item{n}{Range of iterator (starting at \code{1}) for survival functions. This can be the product of number of observations, and number of MCMC samplers in a Bayesian survival function.}
  \item{tout}{Vector of time points for which cumulative incidence functions are requested.}
  \item{Nmax}{Maximum number of subdivisions in the interval [\code{0}, \code{max(tout)}] to be created in the adaptive quadrature algorithm.}
  \item{rel.tol}{Threshold for relative integration error, used as stoppage criterion. It is calculated as the maximum relative error at time point \code{max(tout)} across all causes. Each relative error number is the difference between the Simpson-based and trapezoidal-based numbers, divided by the Simpson-based number.}
  \item{ncores}{Number of parallel threads to use. This is currrently only implemented in \code{C++} mode.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class \code{cfc}, which is a list with the following elements:
\item{ci}{Array of dimensions \code{(length(tout), length(f.list), n)}, cumulative incidence functions for all causes and all values of the iterator, evaluated at all time points indicated by \code{tout}.}
\item{s}{Array of same dimensions as \code{ci}, containing the (unadjusted) survival functions for all causes and all values of the iterator, evaluated at all time points indicated by \code{tout}.}
\item{is.maxiter}{Binary Array of length \code{n}, where \code{1} indicates that subdivision process for quadrature problem applied to survival functions at iteration \code{n} reached maximum set by \code{Nmax} before converging, and \code{0} otherwise.}
\item{n.maxiter}{Number of iterations that did not converge, i.e., \code{sum(is.maxiter)}.}
}
\references{
Haller, B., Schmidt, G., & Ulm, K. (2013). Applying competing risks regression models: an overview. Lifetime data analysis, 1-26.

Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09

Prentice et al (1978). The analysis of failure times in the presence of competing risks. Biometrics, 541-554.
}
\author{
Mansour T.A. Sharabiani, Alireza S. Mahani
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cfc.survreg}}
}
\examples{
\dontrun{

library("survival") # used for constructing survival formulas
library("BSGW") # used for Bayesian survival regression

data("bmt")
# splitting data into training and prediction sets
idx.train <- sample(1:nrow(bmt), size = 0.7 * nrow(bmt))
idx.pred <- setdiff(1:nrow(bmt), idx.train)
nobs.train <- length(idx.train)
nobs.pred <- length(idx.pred)

# prepare data and formula for Bayesian cause-specific survival regression
# using R package BSGW
out.prep <- cfc.prepdata(Surv(time, cause) ~ platelet + age + tcell, bmt)
f1 <- out.prep$formula.list[[1]]
f2 <- out.prep$formula.list[[2]]
dat <- out.prep$dat
tmax <- out.prep$tmax

# estimating cause-specific models
# set nsmp to larger number in real-world applications
nsmp <- 10
reg1 <- bsgw(f1, dat[idx.train, ], control = bsgw.control(iter = nsmp)
  , ordweib = T, print.level = 0)
reg2 <- bsgw(f2, dat[idx.train, ], control = bsgw.control(iter = nsmp)
  , ordweib = T, print.level = 0)

# defining survival function for this model
survfunc <- function(t, args, n) {
  nobs <- args$nobs; natt <- args$natt; nsmp <- args$nsmp
  alpha <- args$alpha; beta <- args$beta; X <- args$X
  idx.smp <- floor((n - 1) / nobs) + 1
  idx.obs <- n - (idx.smp - 1) * nobs
  return (exp(- t ^ alpha[idx.smp] * 
                exp(sum(X[idx.obs, ] * beta[idx.smp, ]))));
}

# preparing function and argument lists
X.pred <- as.matrix(cbind(1, bmt[idx.pred, c("platelet", "age", "tcell")]))
arg.1 <- list(nobs = nobs.pred, natt = 4, nsmp = nsmp
  , alpha = exp(reg1$smp$betas), beta = reg1$smp$beta, X = X.pred)
arg.2 <- list(nobs = nobs.pred, natt = 4, nsmp = nsmp
  , alpha = exp(reg2$smp$betas), beta = reg2$smp$beta, X = X.pred)
arg.list <- list(arg.1, arg.2)
f.list <- list(survfunc, survfunc)

# cause-specific competing-risk
# set rel.tol to smaller number in real-world applications
tout <- seq(from = 0.0, to = tmax, length.out = 10)
out.cfc <- cfc(f.list, arg.list, nobs.pred * nsmp, tout, rel.tol = 1e-2)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
