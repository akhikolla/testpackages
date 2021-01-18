% file stosim/man/ChargeTrain.Rd
% copyright (C) 2014 OpenReliability.org
%
\name{ChargeTrain}
\alias{ChargeTrain}
\title{
  Ficticious system reliability data for an example problem contained in
  "Using Stochastic RAM Analysis to Establish an Optimal Operating Policy".
  Submitted for the "IMA Journal of Management Mathematics, MIMAR 2014 special issue"
}
\description{
  The \code{ChargeTrain} data frame has 3 rows and 11 columns.
}
\usage{
ChargeTrain
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{OpLine}}{
      OpLine: a factor to identify the operating line model applicable to this input.
    }
    \item{\code{EventID}}{
      EventID: a factor to identify the subsystem or EventElement applicable to this input.
    }
    \item{\code{FD}}{
      FD: a factor identifying the distribution for failure data. Only \code{"E"}, \code{"W"},
      and \code{"N"} have been implemented in intended recieving function SimHistory.
    }
    \item{\code{FP1}}{
      A first parameter for the failure distribution (FD). Mean for exponential or normal, 
	  characteristic life for Weibull. 
    }
    \item{\code{FP2}}{
      A second parameter for the failure distribution (FD). Standard deviation for normal, 
	  shape for Weibull.
    }
    \item{\code{FP3}}{
      A third parameter for the failure distribution (FD). A translation parameter applicable to
	  Weibull, only.
    }
    \item{\code{RD}}{
      RD: a factor identifying the distribution for repair time data. Only \code{"L"}, \code{"W"},
      and \code{"N"} have been implemented in intended recieving function SimHistory.
    }
    \item{\code{RP1}}{
      A first parameter for the repair distribution (RD). Meanlog for lognormal, 
	  characteristic life for Weibull, mean for normal. 
    }
    \item{\code{RP2}}{
      A second parameter for the repair distribution (RD). SDlog for lognormal, 
	  shape for Weibull, standard deviation for normal.
    }
    \item{\code{RP3}}{
      A third parameter for the repair distribution (RD). A translation parameter applicable to
	  lognormal or Weibull, only.
    }
    \item{\code{Seed}}{
      An integer seed value for the RNG.
    }
  }
}
\details{
  Typical processing system values were chosen to prepare a likely model utilizing lognormal
  normal and Weibull distributions for repair times (maintainability) as a demonstration.
}
\keyword{datasets}
