\name{calcPredictorOK}
\alias{calcPredictorOK}
\title{
Generate smoothing predictor given smoothing parameters
}
\description{
Assuming that \code{\link{calcGCV}} has been first run to estimate smoothing parameter, this produces a \dQuote{Kriging} predictor of the response.
}
\usage{
calcPredictorOK(FONKgpointls, minKrigPtNbr = blackbox.getOption("minKrigPtNbr"),
                krigmax = NULL, topmode = FALSE, rawPlots = TRUE, cleanResu = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FONKgpointls}{
  Input data frame as produced by \code{\link{buildFONKgpointls}}
}
  \item{minKrigPtNbr}{
  NULL or numeric. At least this many rows (if available) should be selected for Kriging. The default value depends on the number p of predictor variables and is
  90, 159,  500, 1307, 3050, 6560 for p from 1 to 6 (beyond which it is strongly advised to use a non-default value).
}
  \item{krigmax}{
  NULL or Numeric. For large data sets the selected points are not \dQuote{Kriged} all together. Rather, overlapping blocks of rows are selected and are Kriged separately. This sets the size of the blocks. Default depends on the operating system (see source code).
}
  \item{topmode}{
 Controls the way rows are selected. For development purposes, should not be modified
}
  \item{rawPlots}{
  Boolean. Whether to plot one-dimensional \dQuote{profiles} of the raw data.
}
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If \code{""} (the default), print to the standard output connection.
}
}
\value{
Returns invisibly a list with many undocumented elements. Thislist is also stored as a global option \code{"fitobject"}.
}
