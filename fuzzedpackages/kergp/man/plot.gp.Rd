\name{plot.gp}
\alias{plot.gp}

\title{

  Diagnostic Plot for the Validation of a \code{gp} Object

}

\description{

  Three plots are currently available, based on the \code{influence}
  results: one plot of fitted values against response values, one plot
  of standardized residuals, and one qqplot of standardized residuals.

}
\usage{
\S3method{plot}{gp}(x, y, kriging.type = "UK",
    trend.reestim = TRUE, which = 1:3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    
    An object with S3 class \code{"gp"}.

  }
  \item{y}{

    Not used.

  }
  \item{kriging.type}{

    Optional character string corresponding to the GP "kriging" family,
    to be chosen between simple kriging (\code{"SK"}) or universal
    kriging (\code{"UK"}).

  }
  \item{trend.reestim}{

    Should the trend be re-estimated when removing an observation?
    Default to \code{TRUE}.

  }
  \item{which}{

    A subset of {1, 2, 3} indicating which figures to plot (see
    \code{Description} above). Default is 1:3 (all figures).

  }\item{...}{

    No other argument for this method.

  }
}

\details{
  
  The standardized residuals are defined by \eqn{[y(\mathbf{x}_i) -
  \widehat{y}_{-i}(\mathbf{x}_i)] /
  \widehat{\sigma}_{-i}(\mathbf{x}_i)}{(y(xi) - yhat_{-i}(xi)) /
  sigmahat_{-i}(xi)}, where \eqn{y(\mathbf{x}_i)}{y(xi)} is the response at the
  location \eqn{\mathbf{x}_i}{xi},
  \eqn{\widehat{y}_{-i}(\mathbf{x}_i)}{yhat_{-i}(xi)} is the fitted
  value when the \eqn{i}-th observation is omitted (see
  \code{\link{influence.gp}}), and
  \eqn{\widehat{\sigma}_{-i}(\mathbf{x}_i)}{sigmahat_{-i}(xi)} is the
  corresponding kriging standard deviation.

}

\section{Warning}{
  
  Only trend parameters are re-estimated when removing one
  observation. When the number \eqn{n} of observations is small,
  re-estimated values can substantially differ from those obtained with
  the whole learning set.
  
}

\value{

  A list composed of the following elements where \emph{n} is the total
  number of observations.

  \item{mean }{

    A vector of length \emph{n}. The \eqn{i}-th element is the kriging
    mean (including the trend) at the \eqn{i}-th observation number when
    removing it from the learning set.

  }
  \item{sd }{

    A vector of length \emph{n}. The \eqn{i}-th element is the kriging
    standard deviation at the \eqn{i}-th observation number when removing it
    from the learning set.

  }
}


\references{ 

F. Bachoc (2013), "Cross Validation and Maximum Likelihood estimations of
hyper-parameters of Gaussian processes with model
misspecification". \emph{Computational Statistics and Data Analysis},
\bold{66}, 55-69.

N.A.C. Cressie (1993), \emph{Statistics for spatial data}. Wiley series
in probability and mathematical statistics.

O. Dubrule (1983), "Cross validation of Kriging in a unique
neighborhood". \emph{Mathematical Geology}, \bold{15}, 687-699.
  
J.D. Martin and T.W. Simpson (2005), "Use of kriging models to
approximate deterministic computer models". \emph{AIAA Journal},
\bold{43} no. 4, 853-863.
  
M. Schonlau (1997), \emph{Computer experiments and global optimization}.
Ph.D. thesis, University of Waterloo.
  
}

%\author{ O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne. }

\seealso{

  \code{\link{predict.gp}} and \code{\link{influence.gp}}, the
  \code{predict} and \code{influence} methods for \code{"gp"}.

}

\keyword{models}
\keyword{methods}
