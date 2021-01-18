\name{cfc.survreg.survprob}
\alias{cfc.survreg.survprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival probability function for \code{survreg} models
}
\description{
Function for predicting survival probability as a function of time for \code{survreg} regression objects in \pkg{survival} package. It can be used to mix \code{survreg} models with other survival models in competing-risk analysis, using \pkg{CFC} package. This function is used inside \code{cfc.survreg}.
}
\usage{
cfc.survreg.survprob(t, args, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{Time from index. Must be non-negative, but can be a vector.}
  \item{args}{Regression object that is returned by \code{survreg}. If using \code{newdata} for prediction, the \code{x} field of this object must be updated accordingly.}
  \item{n}{Observation index, must be between \code{1} and \code{nrow(args$x)}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Vector of survival probabilities at time(s) \code{t}.
}
\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09
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
library("CFC") # for cfc
data(bmt)
library("randomForestSRC") # for rfsrc
library("survival") # for survreg

prep <- cfc.prepdata(Surv(time, cause) ~ platelet + age + tcell, bmt)
f1 <- prep$formula.list[[1]]
f2 <- prep$formula.list[[2]]
dat <- prep$dat
tmax <- prep$tmax

# building a parametric Weibull regression model
# for cause 1
reg1 <- survreg(f1, dat, x = TRUE) # must keep x for prediction

# building a random forest survival model for cause 2
reg2 <- rfsrc(f2, dat)
# implementing a continuous interface for the random forest
# survival function
rfsrc.survfunc <- function(t, args, n) {
  which.zero <- which(t < .Machine$double.eps)
  ret <- approx(args$time.interest, args$survival[n, ], t, rule = 2)$y
  ret[which.zero] <- 1.0
  return (ret)
}

# constructing function and argument list
f.list <- list(cfc.survreg.survprob, rfsrc.survfunc)
arg.list <- list(reg1, reg2)

# competing-risk analysis
tout <- seq(0.0, tmax, length.out = 10)
# increase rel.tol for higher accuracy
cfc.out <- cfc(f.list, arg.list, nrow(bmt), tout, rel.tol = 1e-3)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
