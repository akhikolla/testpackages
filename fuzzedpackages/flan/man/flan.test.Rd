\name{flan.test}
\alias{flan.test}
\encoding{UTF-8}

\title{Fluctuation Analysis parametric testing}

\description{
  Performs one-sample and two-sample Fluctuation Analysis tests on mutant counts. Returns confidence intervals and p-values, based on asymptotic normality, from the values returned by \code{mutestim}.
}

\usage{

  flan.test(mc, fn = NULL, mfn = NULL, cvfn = NULL,
               fitness = NULL, death = 0., plateff = 1.,
               model = c("LD", "H", "I"), muinf = +Inf,
               mutations0 = 1., mutprob0 = NULL, fitness0 = 1.,
               conf.level = 0.95,
               alternative = c("two.sided", "less", "greater"),
               method = c("ML", "GF", "P0"), winsor = 2000
               )
}

\arguments{
  \item{mc}{ a numeric vector of mutant counts or a list of two numeric vectors of mutant counts.  }
  \item{fn}{ an optional numeric vector of final numbers or a list of two numeric vectors of final numbers. }
  \item{mfn}{ mean final number of cells. Ignored if \code{fn} is non-missing. }
  \item{cvfn}{ coefficient of variation of final number of cells. Ignored if \code{fn} is non-missing. }
  \item{fitness}{ fitness parameter: ratio of growth rates of normal and mutant cells. If \code{fitness} is \code{NULL} (default) then the fitness will be estimated. Otherwise, the given value will
be used to estimate the mean mutation number \code{mutations} }
  \item{death}{ death probability. Must be smaller than 0.5. By default 0. }
  \item{plateff}{ plating efficiency parameter. Must be non-larger than 1. By default 1. Available for \code{GF} method, and for \code{ML} method under \code{LD} model.}
  \item{model}{ statistical lifetime model. Must be one of "LD" (default) for Luria-Delbrück model (exponential lifetimes), "H" for Haldane model (constant lifetimes), or "I" for Inhomogeneous model }
  \item{muinf}{ parameter used only if \code{model} is \code{I}. See details. }
  \item{mutations0}{ null hypothesis value for parameter mutations. See details. }
  \item{mutprob0}{ null hypothesis value for parameter mutprob. See details. }
  \item{fitness0}{ null hypothesis value for parameter fitness. See details. }
  \item{conf.level}{ confidence level of the interval. }
  \item{alternative}{ a character string specifying the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}. }
  \item{method}{ estimation method as a character string: one of \code{ML} (default), \code{P0}, or \code{GF}. See \code{mutestim}. }
  \item{winsor}{ winsorization parameter: positive integer. Only used when \code{method} is \code{ML} or when    \code{method} is \code{P0} and \code{fitness} is NULL. See \code{mutestim}. }
}

\details{
  \code{flan.test} tests the values of parameters \code{mutations}, \code{mutprob}, or \code{fitness}, with \code{mutations0} \code{mutprob0}, \code{fitness0} as null hypotheses.

  If \code{fn} or \code{mfn} are given, mutprob is tested, otherwise, \code{mutations} is tested. If \code{fitness} is given, it is not tested.

  \code{muinf} corresponds to the cumulative division rate on the interval [0 ; +Inf).
  If \code{model} is \code{I}, \code{muinf} has to be finite, else \code{model} is set to "LD".

  \code{alternative} may be a two dimensional vector specifying the alternatives for the two parameters to be tested.

  If \code{mc} is a list, a two-sample test is performed. In that case, the values in \code{mutations0}, \code{mutprob0} and \code{fitness0} apply to the difference between the two samples (with 0 as default value)
  For the two-sample tests, the parameters \code{mfn}, \code{cvfn}, \code{fitness}, \code{death} and \code{plateff} may be two dimensional vectors, where the first (resp. second) components are related to the first (resp. second) sample.
}

\value{
  Returns a list with class \code{"flantest"}. The structure of a \code{"flantest"} object, is similar to that of a \code{"htest"} object (see also \code{\link{t.test}}).
  The class \code{"flantest"} contains the following components :
    \item{Tstat}{ the value of the computed statistic(s). }
    \item{parameter}{ the values of \code{fitness} (if not tested), \code{death} and \code{plateff}. }
    \item{p.value}{ the p-value(s) of the test. }
    \item{conf.int}{ confidence interval(s) for the parameter(s) relative to the specified alternative. }
    \item{estimates}{ the estimate(s). }
    \item{null.value}{ the specified hypothesized value(s). }
    \item{alternative}{ a (vector of) character string(s) describing the alternative hypothesis. }
    \item{model}{ the statistical lifetime model. }
    \item{method}{ method used to compute the estimate(s). }
    \item{data.name}{ a character string giving the name of the complete data. }
}


\seealso{
  \code{\link{mutestim}}.
}

\examples{
# one sample test, without final numbers
b <- unlist(boeal)
# is the mean mutation number greater than 0.7, and the fitness greater than 0.8?
flan.test(b, alternative = "greater", mutations0 = 0.7, fitness0 = 0.8)
# is the mean mutation number less than 0.8, given the fitness?
flan.test(b, alternative = "less", mutations0 = 0.8, fitness = 0.84)

# one sample test, with final numbers
d <- david[[11]]
flan.test(d$mc, d$fn, alternative = "less", mutprob0 = 2e-10, fitness0 = 2)

# two-sample test: test equality of parameters
b1 <- unlist(boeal[1:10])
b2 <- unlist(boeal[11:20])
flan.test(list(b1, b2))

# realistic random sample of size 100: mutation probability 1e-9,
# mean final number 1e9, coefficient of variation on final numbers 0.3,
# fitness 0.9, lognormal lifetimes, 5% mutant deaths
\donttest{
x <- rflan(100, mutprob = 1e-9, mfn = 1e9, cvfn = 0.3, fitness = 0.9, death = 0.05)

# test on mutations and fitness, without final numbers
flan.test(x$mc, mutations0 = 1, fitness0 = 0.9)

# test on mutprob and fitness, with final numbers
flan.test(x$mc, x$fn, mutprob0 = 1e-9, fitness0 = 0.9)

# given fitness

flan.test(x$mc, x$fn, fitness = 0.9, mutprob0 = 1e-9)

# take deaths into account
flan.test(x$mc, x$fn, mutprob0 = 1e-9, fitness0 = 0.9, death = 0.05)

# change method
flan.test(x$mc, x$fn, mutprob0 = 1e-9, fitness0 = 0.9, death = 0.05, method = "GF")

flan.test(x$mc, x$fn, mutprob0 = 1e-9, fitness0 = 0.9, death = 0.05, method = "P0")
# change model
flan.test(x$mc, x$fn, mutprob0 = 1e-9, fitness0 = 0.9, model = "H")

# Two-sample test
y <- rflan(100, mutprob = 1e-9, mfn = 1e9, cvfn = 0.3, fitness = 1.2, death = 0.05)
MC <- list(x$mc, y$mc)
FN <- list(x$fn, y$fn)
#
flan.test(mc = MC, fn = FN, fitness = c(0.9, 1.1), death = 0.05)

}
}
