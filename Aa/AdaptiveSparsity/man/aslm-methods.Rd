\name{aslm-methods}
\alias{logLik.aslm}
\alias{nobs.aslm}
\alias{predict.aslm}
\title{
   Methods handled by lm
}
\description{
These methods are implemented by the \code{\link{lm}} parent class:
\itemize{
  \item{\code{logLik}} { -- Extract log-likelihood }
  \item{\code{predict}} { -- Predict values based on linear model }
  \item{\code{nobs}} { -- Extract the number of observations from a fit }
}
}
%%\usage{
%%\method{logLik}{aslm}(object, ...)
%%\method{predict}{aslm}(object, newdata=NULL, ...)
%%}
%- maybe also 'usage' for other objects documented here.
%%\arguments{
%%  \item{object}{
%%An object of class inheriting from \dQuote{\code{lm}}
%%}
%%  \item{newdata}{
%%An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used
%%}
%%  \item{\dots}{
%%Further arguments.
%%}
%%}
\seealso{
\code{\link{predict.lm}}, \code{\link{logLik.lm}}, \code{\link{nobs}}
}
\keyword{ methods }
\keyword{ models }