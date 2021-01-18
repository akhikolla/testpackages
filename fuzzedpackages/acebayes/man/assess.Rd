\name{assess}
\alias{assess}
\alias{assess.ace}
\alias{assess.pace}

\title{
Compares two designs under the approximate expected utility
}
\encoding{UTF-8}
\description{
Calculates approximations to the expected utility for two designs.
}
\usage{

assess(d1, d2, ...)

\method{assess}{ace}(d1, d2, B = NULL, n.assess = 20, relative = TRUE, ...)

\method{assess}{pace}(d1, d2, B = NULL, n.assess = 20, relative = TRUE, ...)

}

\arguments{
\item{d1, d2}{
\code{d1} should be an object of class \code{"ace"} or \code{"pace"} and \code{d2} should be an object of class \code{"ace"}, \code{"pace"} or \code{"matrix"}. 
}
\item{B}{
An optional argument for controlling the approximation to the expected utility (see \code{\link{ace}}, \code{\link{aceglm}} and \code{\link{acenlm}}). If left unspecified, the value is inherited from the argument \code{d1}.
}
\item{n.assess}{
If \code{d1} was generated from a call to \code{(p)ace} with argument \code{deterministic = FALSE} or from a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{method} being \code{"MC"}, then \code{n.assess} is an optional argument giving the number of evaluations of the stochastic approximation to the expected utility. 
}
\item{relative}{
An optional argument, for when \code{d1} was generated as a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{criterion} being \code{"A"}, \code{"D"} or \code{"E"}, controlling whether the measure of relative efficiency is calculated for \code{d1} relative to \code{d2} (\code{TRUE}; the default) or for \code{d2} relative to \code{d1} (\code{FALSE}).
}
\item{\dots}{
Arguments to be passed to and from other methods.
}
}
\details{
In the case of when \code{d1} was generated from a call to \code{(p)ace} with argument \code{deterministic = FALSE} or from a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{method} being \code{"MC"}, \code{n.assess} evaluations of the stochastic approximation to the expected utility will be calculated for each of the designs from \code{d1} and \code{d2}. Otherwise, one evaluation of the deterministic approximation to the expected utility will be calculated for each of the designs from \code{d1} and \code{d2}.

In the case when \code{d1} was generated as a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{criterion} being \code{"A"}, \code{"D"} or \code{"E"}, the relative D-, E-, or A-efficiency of the two designs will be calculated. The direction of the relative efficiency can be controlled by the \code{relative} argument.
}

\value{
The function will an object of class \code{"assess"} which is a list with the following components:

\item{U1}{
In the case of when \code{d1} was generated from a call to \code{(p)ace} with argument \code{deterministic = FALSE} or from a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{method} being \code{"MC"}, \code{U1} will be a vector of \code{n.assess} evaluations of the stochastic approximation to the expected utility for design \code{d1}. Otherwise, \code{U1} will be a scalar of one evaluation of the deterministic approximation to the expected utility for design \code{d1}.
}
\item{U2}{
In the case of when \code{d1} was generated from a call to \code{(p)ace} with argument \code{deterministic = FALSE} or from a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{method} being \code{"MC"}, \code{U2} will be a vector of \code{n.assess} evaluations of the stochastic approximation to the expected utility for design \code{d2}. Otherwise, \code{U2} will be a scalar of one evaluation of the deterministic approximation to the expected utility for design \code{d2}.
}
\item{eff}{
In the case when \code{d1} was generated as a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{criterion} being \code{"A"}, \code{"D"} or \code{"E"}, \code{eff} is a scalar of the relative D-, E-, or A-efficiency of the two designs. Otherwise it will be \code{NULL}. 
}
\item{d1}{
The argument \code{d1}.
}
\item{d2}{
The argument \code{d2}.
}
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\seealso{
\code{\link{ace}}, \code{\link{pace}}, \code{\link{aceglm}}, \code{\link{acenlm}}, \code{\link{paceglm}}, \code{\link{pacenlm}}.
}
\examples{
## This example involves finding a Bayesian D-optimal design for a 
## compartmental model with n = 18 runs. There are three parameters. 
## Two parameters have uniform priors and the third has a prior 
## point mass. 

n <- 18
k <- 1
p <- 3
set.seed(1)
start.d <- randomLHS(n = n, k = k) * 24
colnames(start.d) <- c("t")

a1<-c(0.01884, 0.298)
a2<-c(0.09884, 8.298)

prior <- list(support = cbind(rbind(a1, a2), c(21.8, 21.8)))
colnames(prior[[1]]) <- c("theta1", "theta2", "theta3") 

example <- acenlm(formula = ~ theta3 * (exp( - theta1 * t) - exp( - theta2 * t)), 
start.d = start.d, prior = prior, lower = 0, upper = 24, N1 = 2, N2 = 0)

## Compute efficiency of final design compared to starting design.
assess(d1 = example, d2 = start.d)

## Should get 

# Approximate expected utility of d1 = 15.40583 
# Approximate expected utility of d2 = 11.26968 
# Approximate relative D-efficiency = 396.9804% 

}
