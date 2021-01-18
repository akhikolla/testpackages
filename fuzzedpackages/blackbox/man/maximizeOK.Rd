\name{maximizeOK}
\alias{maximizeOK}
\title{
Find maximum of predicted response surface
}
\description{
Assuming that \code{\link{calcPredictorOK}} has been first run to produce a predictor of  the response surface, this finds its constrained maximum in the convex envelope of the smoothed data.
}
\usage{
maximizeOK(fitobject = blackbox.getOption("fitobject"), cleanResu = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitobject}{
  Return object of \code{\link{calcPredictorOK}}.
}
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If \code{""} (the default), print to the standard output connection.
}
}
\value{
A list with element
\item{par }{predictor values maximizing the predicted response (in the parameter space used for Kriging)}
\item{value }{maximum of the predicted response}
\item{canonVP}{Representation of \code{par} in canonical parameter space}
and possibly other elements (i) returned by an optimization function such as \code{optim}; (ii) values of additional transformed parameters; (iii) cryptic information whether maximization occurred at some boundary of the convex envelope.
%% ...
}
