\name{calc1Dprofiles}
\alias{calc1Dprofiles}
\alias{calcProfileLR}
\alias{calc2D3Dplots}
\title{
One and two-dimensional profiles, and surface plots
}
\description{
Assuming that \code{\link{calcPredictorOK}} and \code{\link{maximizeOK}} have been first run:
 \code{calc1Dprofiles} plots 1D profiles of a predicted likelihood surface for each of the parameters. Poor profiles mayresult when only local optima are found for some parameter values. The next function provides an improvement over this.
 \code{calcProfileLR} plots 2D profiles of the predicted response surface relative to its maximum for pairs of parameters. It also prots 1D profiles taking benefit of the computation effort for the 2D profiles.
 \code{calc2D3Dplots} plots the predicted response surface (no profile) in different ways depending on the number of parameters.

 These functions have almost no arguments, as almost all control is through global controls. See in particular \code{gridStepsNbr} (for profile plots) and \code{graphicPars} in \code{\link{blackbox.options}}.
}
\usage{
calc1Dprofiles(varNames=blackbox.getOption("spec1DProfiles"))
calcProfileLR(varNames=blackbox.getOption("fittedNames"),
              pairlist=list(),
              cleanResu="")
calc2D3Dplots(plotFile=NULL,pairlist=list())
}
\arguments{
  \item{plotFile}{If a character string, the name of the file where plots are written. Otherwise, plots are output to the screen.}
  \item{varNames}{A character vector specifying the names of predictor variables to be considered. For \code{calc1Dprofiles} (used in conjunction with the \code{Migraine} software), if the default argument is \code{NULL}, all variable canonical parameters plus some composite ones may be considered (see the source code for details).}
  \item{pairlist}{A list of character vectors. Each vector describes a pair of predictor variables. With the default value \code{list()}, a default non-empty list may be constructed when \code{calc2D3Dplots} or \code{calcProfileLR} is typically used in conjunction with the \code{Migraine} software (see the source code for details).}
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If it is \code{""} (the default), print to the standard output connection.
}
}
\details{
If there is only \bold{one} parameter, \code{calc2D3Dplots} plots the predicted response as function of this parameter

If there are \bold{two} parameters, \code{calc2D3Dplots} plots the response surface both as a 2D surface plot and as a 3D perspective plot, and \code{calcProfileLR} also produces a plot of the response surface (no profiling is needed) relative to its maximum (hence, a likelihood ratio, if the response is a likelihood).

If there are \bold{more} parameters, \code{calc2D3Dplots} plots a \dQuote{slice} of the predicted surface, both as a 2D surface plot and as a 3D perspective plot, for each pair of parameters. A slice plot for a pair of parameters fixes all other parameters to values maximizing   the response (hence, maximum likelihood estimates, if the response is a likelihood). \code{calcProfileLR} plots the profile response surface relative to its maximum (hence, a profile likelihood ratio, if the response is a likelihood) for pairs of parameters in \code{varNames}.

Two dimensional profile plots not only require many numerical maximizations, but will look ugly whenever one of these maximizations fails to find the right maximum, hence additional intensive computations are performed to minimize this problem. As a result, they are quite slow to compute, unless a low \code{gridStepsNbr} (say < 16) is used, in which case they do not look smooth.
}
\value{
Returns NULL invisibly
}
