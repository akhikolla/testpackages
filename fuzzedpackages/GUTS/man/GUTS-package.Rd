\encoding{UTF-8}

\name{GUTS-package}
\alias{GUTS-package}
%\alias{attr<-}
\alias{attr<-.GUTS}
%\alias{attributes<-}
\alias{attributes<-.GUTS}
%\alias{mostattributes<-}
\alias{mostattributes<-.GUTS}
\alias{print.GUTS}
\alias{modguts}
\alias{[[<-.GUTS}
\alias{$<-.GUTS}

\docType{package}

\title{Fast Calculation of the Likelihood of a Stochastic Survival Model}

\description{GUTS (General Unified Threshold model of Survival) is a stochastic survival model for ecotoxicology.  The package allows for the definition of exposure and survival time series as well as parameter values, and the fast calculation of the survival probabilities as well as the logarithm of the corresponding likelihood.}

\details{
\tabular{ll}{
Package: \tab GUTS\cr
Type: \tab Package\cr
%Version: \tab 1.1.0\cr
%Date: \tab 2019-01-29\cr
License: \tab GPL (>= 2)\cr
}
A GUTS object is a special list of class \dQuote{GUTS}.  Functions \code{\link{guts_setup}}, \code{\link{guts_calc_loglikelihood}} and \code{\link{guts_calc_survivalprobs}} are available to create and work with GUTS objects.  A data set \link{diazinon} is also included.  See links for more details.
}

\author{Carlo Albert \email{carlo.albert@eawag.ch}, Sören Vogel \email{soeren.vogel@posteo.ch}, Oliver Jakoby \email{oliver.jakoby@rifcon.de}, Alexander Singer \email{alexander.singer@rifcon.de} and Dirk Nickisch \email{dirk.nickisch@rifcon.de}

Maintainer: Oliver Jakoby \email{oliver.jakoby@rifcon.de}}

\references{Albert, C., Vogel, S., and Ashauer, R. (2016). Computationally efficient implementation of a novel algorithm for the General Unified Threshold Model of Survival (GUTS). PLOS Computational Biology, 12(6), e1004978. doi: 10.1371/journal.pcbi.1004978.

Jager, T., Albert, C., Preuss T., and Ashauer R. (2011). General Unified Threshold Model of Survival -- a toxicokinetic toxicodynamic framework for ecotoxicology. Environmental Science \& Technology, 45(7), 2529--2540, doi: 10.1021/es103092a.}

\seealso{\code{\link{guts_setup}}, \code{\link{guts_calc_loglikelihood}}, \code{\link{guts_calc_survivalprobs}}, \code{\link{guts_report_damage}}, \code{\link{diazinon}}, \code{\link{Rcpp}}}
