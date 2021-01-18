\name{number of mutant cells distribution}
\alias{dflan}
\alias{pflan}
\alias{qflan}
\alias{rflan}
\encoding{UTF-8}


\title{The distribution of mutant cell counts}

\description{
  Density, distribution function, quantile function and random generation for mutant cell counts.
}

\usage{
  dflan(m, mutations = 1., fitness = 1., death = 0., plateff = 1.,
        model = c("LD", "H", "I"), muinf = +Inf)
  pflan(m, mutations = 1., fitness = 1., death = 0., plateff = 1.,
        model = c("LD", "H", "I"), muinf = +Inf, lower.tail = TRUE)
  qflan(p, mutations = 1., fitness = 1., death = 0., plateff = 1.,
        model = c("LD", "H", "I"), muinf = +Inf, lower.tail = TRUE)
  rflan(n, mutations = 1., mutprob = NULL, fitness = 1., death = 0., plateff = 1.,
        dist = list(name = "lnorm", meanlog = -0.3795851, sdlog = 0.3016223),
        distfn = NULL,
        mfn = 1e9, cvfn = 0,
        muih = list(mu = NULL, muinv0 = NULL),
        ...)
}

\arguments{
  \item{m}{ vector of quantiles. }
  \item{p}{ vector of probabilities. }
  \item{n}{ number of observations. If \code{length(n) > 1}, the length is taken to be the number required. }
  \item{mutations}{mean number of mutations: positive integer. Ignored if \code{mutprob} is non-empty. }
  \item{mutprob}{ mutation probability: numeric between 0 and 1. By default empty. See details. }
  \item{fitness}{ fitness parameter: numeric positive. By default 1. }
  \item{death}{ death probability: numeric between 0 and 0.5. By default 0. }
  \item{plateff}{ plating efficiency parameter: numeric number between 0 and 1. By default 1. See Details. }
  \item{dist}{ lifetime distribution for mutant cells. See Details. }
  \item{model}{ statistical lifetime model. Must be one of "LD" (default) for Luria-Delbrück model (exponential lifetimes), "H" for Haldane model (constant lifetimes), or "I" for Inhomogeneous model }
  \item{distfn}{ final number of cells distribution. Must be one of "lnorm" (default if \code{cvfn} is positive) or "gamma". }
  \item{mfn}{ mean final number of cells: numeric positive. }
  \item{cvfn}{ coefficient of variation of final numbers of cells: numeric, default 0. If non-zero and if \code{mutprob} is empty, compute the sample with \code{mutprob} = \code{mutations}/\code{mfn}}
  \item{muih}{ functions for inhomogeneous model. See details. }
  \item{muinf}{ parameter used only if \code{model} is \code{I}. See details. }
  \item{lower.tail}{ logical: if TRUE (default), probabilities are P[X <= x] otherwise, P[X > m]. }
  \item{...}{ Arguments to be passed for \code{muih}, such that hyperparameters. }
}


\details{
  The argument \code{dist} is a list beginning with the distribution name followed by its parameters, and must be one the 4 following distributions: "dirac", "exp", "lnorm"(meanlog, sdlog), "gamma"(shape, scale).
  Note that the parameters related to the "dirac" and "exp" cases are directly computed with inputs \code{fitness} and \code{death}.

  If \code{muih$mu} is non-empty, a time dependent growth model is performed.
  The function \code{muih$mu} is a bivariate function which represents roughly the mean growth of clone durig an time interval.
  The function \code{muih$muinv0} is the inverse function of \code{muih$mu(0,.)}.
  If \code{muih$muinv0} is empty, it is approximated using \code{muih$mu}.

  \code{muinf} corresponds to the cumulative division rate on the interval [0 ; +Inf).
  It is the limit as t tends to infinity of \code{muih$mu(0,t)}.
  If \code{model} is \code{I}, \code{muinf} has to be finite, else \code{model} is set to "LD"

  If \code{mutprob} is non-empty and if \code{cvfn} is zero, the sample of mutants is computed with \code{mutations} as the product of \code{mutprob} by \code{mfn}.
  If \code{cvfn} is positive, the final numbers of cells are generated with the log-normal (if \code{distfn} is "lnorm" or empty) or the gamma (if \code{distfn} is "gamma") distribution with mean \code{mfn} and coefficient of variation \code{cvfn}.
  The sample of mutants is then generated conditionally to the sample of final numbers of cells.

  A plating efficiency \code{plateff} smaller than 1 can be taken into account in the \code{dflan, pflan} and \code{qflan} functions only under \code{LD} model.

  If the plating efficiency \code{plateff} is smaller than 1 in the \code{flan} function, then each mutant count \code{mc} is replaced by a binomial number with parameters \code{size = mc} and \code{prob = plateff}.

}


\value{

    \code{dflan} gives the density, \code{pflan} gives the distribution
     function, \code{qflan} gives the quantile function, and \code{rflan}
     generates a random sample.

     \code{rflan} returns a list with two arguments, each with length \code{n}: a vector of integers
     \code{$mc} (mutant counts), and a vector of numeric \code{$fn} (final numbers of cells).
}


\seealso{
  \code{link{draw.clone}}
}

\examples{
#----------------------- distributions ----------------------------------

# Luria-Delbrück model, mean number of mutations 1, fitness parameter 1
dflan(0:10)
pflan(0:10)
qflan(c(0.95, 0.99))

# Luria-Delbrück model, mean number of mutations 2, fitness parameter 0.5
qflan(c(0.95, 0.99), mutations = 2, fitness = 0.5)
\donttest{
qflan(c(0.05, 0.01), mutations = 2, fitness = 0.5, lower.tail = FALSE)

# Haldane model, mean number of mutations 2, fitness parameter 0.5
qflan(c(0.95, 0.99), mutations = 2, fitness = 0.5, model = "H")
}
#---------------------- random samples ----------------------------------
# lognormal lifetime distribution, lognormal final numbers
X <- rflan(100, cvfn = 0.3)
X$mc              # mutant counts
X$fn              # final numbers

# lognormal lifetime distribution, gamma final numbers
X <- rflan(100, distfn = "gamma", cvfn = 0.3)
X$mc              # mutant counts
X$fn              # final numbers

# mean number of mutations 2, fitness parameter 0.5 (realistic model, but slow)
\donttest{
# lognormal lifetimes
rflan(1000, mutations = 2, fitness = 0.5)$mc
# gamma lifetimes
rflan(1000, mutations = 2, fitness = 0.5, dist = list(name = "gamma", shape = 10, scale = 0.05))$mc
}
# exponential lifetimes (Luria-Delbrück model, fast)
rflan(1000, mutations = 2, fitness = 0.5, dist = list("exp"))$mc

# constant lifetimes (Haldane model, fast)
rflan(1000, mutations = 2, fitness = 0.5, dist = "dirac")$mc

# specifying mutation probability and mean final number
\donttest{
rflan(1000, mutprob = 2e-9, mfn = 2e9, fitness = 0.5, dist = list("dirac"))$mc


# positive cell death probability
rflan(1000, mutprob = 2e-9, mfn = 2e9, death = 0.1, fitness = 0.5, dist = list("exp"))$mc
}
#------------- Inhomogeneous model ---------
f <- function(t,n0,nf) nf/(1+(nf/n0-1)*exp(-t))  # Logistic function: mean growth of a clone
finv <- function(s,n0,nf) -log((nf/s-1)/(nf/n0-1))  # Inverse function of f

mu <- function(s,t,n0,nf) log(f(t,n0,nf)/f(s,n0,nf))   # Definition of mu
muinv0 <- function(u,n0,nf) finv(f(0.,n0,nf)*exp(u),n0,nf)   # Inverse function of mu

rflan(1000, muih = list(mu = mu, muinv0 = muinv0), n0 = 1e2, nf = 1e6)$mc

}
