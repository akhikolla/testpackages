\name{summary.cfc}
\alias{summary.cfc}
\alias{plot.summary.cfc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing and plotting output of \code{cfc}
}
\description{
\code{summary} method for class \code{\link{cfc}}.
}
\usage{
\method{summary}{cfc}(object
  , f.reduce = function(x) x
  , pval = 0.05, ...)
\method{plot}{summary.cfc}(x, which = c(1, 2), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class "cfc", usually the result of a call to \code{\link{cfc}}.}
  \item{f.reduce}{Function to be applied to each sub-array of \code{object$ci} (cumulative incidence) and \code{object$s} (survival probability).}
  \item{pval}{Desired significance level for confidence intervals produced by \code{summary.cfc}. We essentially set the argument \code{probs} to \code{c(pval/2, 0.5, 1-pval/2)} when calling \code{quantile}.}
  \item{x}{An object of class "summary.cfc", usually the result of a call to \code{\link{summary.cfc}}.}
  \item{which}{Vector of integers, indicating which plot(s) must be produced: 1) cumulative incidence functions, one per cause. For each cause, median and credible bands are plotted vs. time-from-index. 2) (unadjusted) survival functions, one per cause. Similar to (1), median and credible bands are plotted.}
  \item{...}{Further arguments to be passed to \code{f.reduce} (for \code{summary.cfc}).}
}
\value{
Recall that the survival probability and cumulative incidence arrays returned by \code{cfc} are three-dimensional, and their first two dimensions indicate 1) time points and 2) causes. \code{f.reduce} is expected to produce an array of a fixed length, when applied to each sub-array, \code{ci[i, j, ]} and \code{s[i, j, ]}. The end-result is two three-dimensional array, where the first two dimensions are identical to its input arrays. This 3D array is then passed to the \code{quantile} function to compute median and credible bands. There is a special case where \code{f.reduce} returns a scalar, rather than an array, when applied to each sub-array. In this case, quantile calculation is meaningless and we return simply these point estimates. In summary, the return object from \code{summary} is a list with elements: 1) \code{ci} (cumulative incidence), 2) \code{s} (survival), and 3) \code{quantiles}, a boolean flag indicating whether the cumulative incidence and survival arrays returned are quantiles or point estimates.
}

\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09
}

\author{Alireza S. Mahani, Mansour T.A. Sharabiani}

\seealso{
\code{\link{cfc}}, \code{\link{summary}}.
}

\examples{
\dontrun{

library("BSGW") # used for Bayesian survival regression

data(bmt)
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
out.cfc <- cfc(f.list, arg.list, nobs.pred * nsmp, tout, rel.tol = 1e-2)

# summarizing (and plotting) the results
# this function calculates the population-average CI and survival, one
# per each MCMC sample; therefore, the quantiles produced by the summary
# method, correspondingly, reflect our confidence in population-average values
my.f.reduce <- function(x, nobs, nsmp) {
  return (colMeans(array(x, dim = c(nobs, nsmp))))
}
my.summ <- summary(out.cfc, f.reduce = my.f.reduce, nobs = nobs.pred, nsmp = nsmp)

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

