\encoding{UTF-8}
% Last updated 2019-07-04

\name{GUTS}

\alias{GUTS}
\alias{guts}
\alias{guts_setup}
\alias{guts_calc_loglikelihood}
\alias{guts_calc_survivalprobs}
\alias{guts_report_damage}
\alias{guts_report_sppe}
\alias{guts_report_squares}



\title{Fast Calculation of the Likelihood of a Stochastic Survival Model}



\description{GUTS (General Unified Threshold model of Survival) is a stochastic survival model for ecotoxicology. The package allows for the definition of exposure and survival time series as well as parameter values, and the fast calculation of the survival probabilities as well as the logarithm of the corresponding likelihood. 

The package implements the GUTS-SIC (also called GUTS-RED) variants that assume a one-compartment model with first-order toxicokinetics.}


\usage{
guts_setup(C, Ct, y, yt, dist = "lognormal",
	model = "Proper", N = 1000, M = 10000)

guts_calc_loglikelihood(gobj, par, external_dist = NULL)

guts_calc_survivalprobs(gobj, par, external_dist = NULL)

guts_report_damage(gobj)

guts_report_sppe(gobj)

guts_report_squares(gobj)
}


\arguments{%
	\item{C}{Numeric vector of concentrations.  Vector must contain at least 2 values and be of the same length as \code{Ct}.%
	}
	\item{Ct}{Numeric vector of concentration time points.  Vector must contain at least 2 values and be of the same length as \code{C}. Time points must start at 0, and contain unique values in ascending order.%
	}
	\item{y}{Integer vector (counts) of survivors.  Vector must contain at least 2 values and be of the same length as \code{yt}.  \code{y} must not be ascending.%
	}
	\item{yt}{Numeric vector of survivor time points.  Vector must contain at least 2 values and be of the same length as \code{y}. Time points must start at 0, and contain unique values in ascending order.  Survivor information at time points later than the latest concentration time point will be disregarded (with a warning).%
	}
	\item{dist}{Distribution as character, either \dQuote{lognormal} (default), \dQuote{loglogistic}, \dQuote{external} or \dQuote{delta}.%
	}
	\item{model}{Model as character, either \dQuote{Proper} (for full model, the default), \dQuote{IT} (for individual tolerance), or \dQuote{SD} (for stochastic death).%
	}
	\item{N}{Integer.  Thresholds sample length.  Must be greater than 2.%
	}
	\item{M}{Integer.  Number of time grid points.  Must be greater than 1.%
	}
	\item{gobj}{GUTS object.  The object to be updated (and used for the calculation).%
	}
	\item{par}{Numeric vector of parameters.  See details below.%
	}
	\item{external_dist}{Numeric vector containing the distribution of individual thresholds. Only used if \code{dist = 'external'}. See details below.%
	}
} % End of \arguments



\details{%


\subsection{Functions}{%
Use \code{guts_setup} to define (or alter) a GUTS object. Various checks are applied to the data. On success, a GUTS object will be created.

Use \code{guts_calc_loglikelihood} to calculate the survival probabilities and the corresponding loglikelihood for a given set of parameters.  The function is very fast and can be used in routines for parameter estimation.  The function returns the loglikelihood, however it also updates the fields \code{par}, \code{S}, \code{D}, \code{SPPE}, \code{squares}, \code{zt} and \code{LL} of the GUTS-object.

\code{guts_calc_survivalprobs} is a convenience wrapper that can be used for predictions; it returns the survival probabilities, however it also updates the fields \code{par}, \code{S}, \code{D}, \code{SPPE}, \code{squares}, \code{zt} and \code{LL} of the GUTS-object.

\code{guts_report_damage} returns a data.frame with time grid points and the damage for each of these. The function reports the damage that was calculated in the previous call to \code{guts_calc_loglikelihood} or \code{guts_calc_survivalprobs}.

\code{guts_report_squares} returns the sum of squares. The function reports the sum of squares that was calculated in the previous call to \code{guts_calc_loglikelihood} or \code{guts_calc_survivalprobs}.

\code{guts_report_sppe} returns the survival-probability prediction error (SPPE). The function reports the SPPE that was calculated in the previous call to \code{guts_calc_loglikelihood} or \code{guts_calc_survivalprobs}.
}


\subsection{ Models, Parameters, and Distributions}{%

The GUTS package provides three model types:
\itemize{%
	\item Proper: a GUTS-SIC-Proper (also called GUTS-RED-Proper) model using random individual tolerances and a stochastic death process, when individual tolerances are exceeded.
	\item IT: a GUTS-SIC-IT (GUTS-RED-IT) individual tolerance model using random individual tolerances. If an individual's tolerance threshold is exceeded, the individual dies.
	\item SD: a GUTS-SIC-SD (GUTS-RED-SD) stochastic death model using a stochastic death process above a population-wide tolerance threshold. The tolerance-threshold is the same for all individuals.
}


The Proper GUTS model requires the following parameters \code{par}, while variants IT and SD are based on a reduced subset (as indicated in brackets). Parameter values in \code{par} must be ordered as listed here:
\itemize{%
	\item hb: background mortality rate (Proper, IT, SD)
	\item ke: dominant rate constant (Proper, IT, SD)
	\item kk: killing rate (Proper, SD)
	\item further parameters for the tolerance threshold (in SD) or the threshold distribution \code{dist} (in Proper and IT) 
}

For model type \dQuote{SD} (stochastic death), required parameters \code{par[1:4]} are \code{kb}, \code{ke}, \code{kk} and \code{mn}, which is the population-wide tolerance threshold. For backwards compatibility this model type can be initiated setting \code{dist = "Delta"} and \code{model = "Proper"}.
For model type \dQuote{IT} (individual tolerance), required parameters \code{par[1:2]} are \code{kb}, \code{ke}, as well as respective distribution parameters (continued from \code{par[3]}). Parameter (\code{kk}) is set internally to infinity and must not be provided.
For model type \dQuote{Proper}, all parameters are needed. \code{par[1:3]} take \code{kb}, \code{ke}, \code{kk}, distribution parameters follow from \code{par[4]}).

For model types \dQuote{Proper} and \dQuote{IT} individual tolerance thresholds are created internally. Individual tolerances are drawn from the specified distribution \code{dist}:
\itemize{%
	\item "lognormal": requires the parameters \code{mn} and \code{sd} which are the mean and standard deviation of the lognormal random distribution. In contrast to parameters \code{meanlog} and \code{sdlog} of function \code{\link{dlnorm}}, these parameters are not on the logscale. They relate in the following way:
\deqn{sdlog = \sqrt{\ln{\frac{1 + sd^2}{mn^2}}}}{%
      sdlog = (ln( 1 + sd^2 / mn^2))^0.5}
\deqn{meanlog = \ln{mn} - \frac{1}{2} * sdlog^2}{%
      meanlog = ln(mn) - 0.5 * sdlog^2}
	\item "loglogistic": requires the parameters \eqn{mn = scale = median} and \eqn{beta = shape}.
	\item "external": uses random variates provided to \code{external_dist}. With this option GUTS can be run with arbitrarily distributed individual tolerance thresholds. With the option \dQuote{external} only parameters \code{hb}, \code{ke} and \code{kk} are required. Further, the thresholds sample length \code{N} is internally adjusted to the length of the external vector of random variates \code{external_dist}. The adjustment of \code{N} is notified by a warning.
}
For performance reasons the implemented distributions \dQuote{lognormal}  and \dQuote{loglogistic} are approximated using importance sampling. The option \dQuote{external} generally performance well, but might require a larger thresholds sample (i.e. \code{length(external_dist)} should be large). 

The number of parameters is checked according to \code{dist} and \code{model}.  Wrong number of parameters invoke an error, wrong parameter values (e.g., negative values) invoke a warning, and the loglikelihood is set to \code{-Inf}.

} % End of \subsection{ Models, Parameters, and Distributions}.


\subsection{Field and Attribute Access}{%
Fields and attributes of an object of class \dQuote{GUTS} are read-only.  To prevent accidental change of fields or attributes, replacement functions were rewritten throwing an error when used.  Always use function \code{guts_setup} to create objects or modify fields on existing objects. Functions \code{guts_calc_loglikelihood} and \code{guts_calc_survivalprobs} update an object's fields \code{par} (parameters), \code{D} (damage), \code{squares} (sum of squares), \code{SPPE} (survival-probability prediction error), \code{S} (survival probabilities) and \code{LL} (the loglikelihood).
} % End of \subsection{Field and Attribute Access}.
} % End of \details



\value{
\code{guts_setup} returns a list of class \dQuote{GUTS} with the following fields:

\item{C}{Concentrations.}
\item{Ct}{Concentration time points.}
\item{y}{Survivors.}
\item{yt}{Survivor time points.}
\item{dist}{Distribution.}
\item{model}{Model.}
\item{N}{Sample length.}
\item{M}{Time grid points.}
\item{par}{Parameters.}
\item{S}{Vector of survivor probabilities.}
\item{D}{Vector of internal damage for each of the \code{M} time grid points.}
\item{squares}{Sum of squares}
\item{SPPE}{Survival-probability prediction error.}
\item{LL}{The loglikelihood.}

%The object has the following attributes (for internal use):
%
%\item{class}{\code{GUTS}}
%\item{experiment}{Internal representation of distribution and model.}
%\item{wpar}{Internal representation of parameters.}
%\item{par_pos}{Internal representation of parameter positions.}

\code{guts_calc_loglikelihood} returns the loglikelihood.

\code{guts_calc_survivalprobs} returns the survival probabilities.

\code{guts_report_damage} returns the damage.

\code{guts_report_squares} returns the sum of squares.

\code{guts_report_sppe} returns the survival-probability prediction error (SPPE).

} % End of \value.



\note{%
The GUTS project web site can be found here: \href{http://guts.r-forge.r-project.org}{http://guts.r-forge.r-project.org}.  For questions and discussion, please subscribe to the mailing list there.
} % End of \note.



\references{Albert, C., Vogel, S., and Ashauer, R. (2016). Computationally efficient implementation of a novel algorithm for the General Unified Threshold Model of Survival (GUTS). PLOS Computational Biology, 12(6), e1004978. doi: 10.1371/journal.pcbi.1004978.

Jager, T., Albert, C., Preuss, T., and Ashauer, R. (2011). General Unified Threshold Model of Survival -- a toxicokinetic toxicodynamic framework for ecotoxicology. Environmental Science \& Technology, 45(7), 2529--2540, doi: 10.1021/es103092a.

Ashauer, R., Albert, C., Augustine, S., Cedergreen, N., Charles, S., Ducrot, V., Focks, A., Gabsi, F., Gergs, A., Goussen, B., Jager, T., Kramer, N.I., Nyman, A.-M., Poulsen, V., Reichenberger, S., Schäfer, R.B., Van den Brink, P.J., Veltman, K., Vogel, S., Zimmer, E.I., Preuss, T.G. (2016) Modelling survival: exposure pattern, species sensitivity and uncertainty. Scientific Reports, 6, 1.

Jager, T., Ashauer, R. (2018). Modelling survival under chemical stress. A comprehensive guide to the GUTS framework. Leanpub: \url{https://leanpub.com/guts_book}, \url{http://www.debtox.info/book_guts.html}

EFSA PPR Panel (EFSA Panel on Plant Protection Products and their Residues), Ockleford, C., Adriaanse, P., Berny, P., Brock, T., Duquesne, S., Grilli, S., Hernandez-Jerez, A.F., Bennekou, S.H., Klein, M., Kuhl, T., Laskowski, R., Machera, K., Pelkonen, O., Pieper, S., Smith, R.H., Stemmer, M., Sundh, I., Tiktak, A., Topping, C.J., Wolterink, G., Cedergreen, N., Charles, S., Focks, A., Reed, M., Arena, M., Ippolito, A., Byers, H. and Teodorovic, I. (2018). Scientific Opinion on the state of the art of Toxicokinetic/Toxicodynamic (TKTD) effect models for regulatory risk assessment of pesticides for aquatic organisms. EFSA Journal, 16(8):5377, 188 pp. \url{https://doi.org/10.2903/j.efsa.2018.5377}
}


\author{Carlo Albert \email{carlo.albert@eawag.ch}, Sören Vogel \email{soeren.vogel@posteo.ch}, Oliver Jakoby \email{oliver.jakoby@rifcon.de}, Alexander Singer \email{alexander.singer@rifcon.de} and Dirk Nickisch \email{dirk.nickisch@rifcon.de}

Maintainer: Oliver Jakoby \email{oliver.jakoby@rifcon.de}}


\seealso{\code{\link{diazinon}}, \code{\link{GUTS-package}} and the package vignettes for examples on how to calibrate and project GUTS-models.}



\examples{
data(diazinon)

# create GUTS object to calculate the Proper model 
#  using a log-normal distribution of tolerance thresholds  
gts.lognormal <- guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "lognormal", model = "Proper")

# calculate likelihood of Proper model using log-normal distribution
guts_calc_loglikelihood(  
  gts.lognormal, 
  c(0.051, 0.126, 1.618, 19.099, 6.495))
gts.lognormal  # show GUTS object

# repeating calculation above 
#  with threshold values from an external log-normal distribution. 
#  Note, we need to account for the different parametrisations
#  used in the GUTS-package and in rlnorm
sigma2 <- log( 1 + 6.495^2 / 19.099^2)
mu <- log(19.099) - 0.5 * sigma2
lognormal.thresholds <- rlnorm(1000, meanlog = mu, sdlog = sqrt(sigma2))
gts.external <- guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "external", model = "Proper")
guts_calc_loglikelihood(
  gts.external, 
  c(0.051, 0.126, 1.618), external_dist = lognormal.thresholds) 
# -> Results using external and internal distributions are comparable

# create GUTS object to calculate the Proper model 
#  using a log-logistic distribution of tolerance thresholds
gts.loglogistic <- guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "loglogistic", model = "Proper")
guts_calc_survivalprobs( # returning survival probabilities
  gts.loglogistic, 
  c(0.01, 0.2, 0.3, 3, 2)) 

str(guts_report_damage(gts.loglogistic)) # returning damage

# calculate survival probabilities with IT model 
#  using a log-logistic distribution of tolerance thresholds
guts_calc_survivalprobs(
  guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "loglogistic", model = "IT"), 
  c(0.01, 0.2, 3, 2)) 

# calculate survival probabilities with an SD model with a fixed tolerance threshold
guts_calc_survivalprobs(
  guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "loglogistic", model = "SD"), 
  c(0.01, 0.2, 0.3, 3))


\dontrun{guts_calc_survivalprobs(gts.external, rep(.5, 3))} 
# Warning and no result, because no external distribution was specified

\dontrun{guts_calc_survivalprobs(gts.loglogistic , 1:4 )} # Error.

\dontrun{gts.loglogistic[["C"]] <- 1:3} # Error.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line



