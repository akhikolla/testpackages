% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edina.R
\name{summary.edina}
\alias{summary.edina}
\title{Summarize the EDINA Object}
\usage{
\method{summary}{edina}(object, alpha = 0.05, ...)
}
\arguments{
\item{object}{An \code{edina} object}

\item{alpha}{Defining region to indicate the level of extremeness
the data must before the model is problematic.}

\item{...}{Not used.}
}
\value{
A summary object that includes everything in the original \code{\link[=edina]{edina()}} object
and:
\itemize{
\item \code{model_fit}: Matrix of model fit summary statistics.
\item \code{alpha}: Alpha-value used to compute \code{\link[=PPP]{PPP()}}s.
}
}
\description{
Provide a more detailed view inside of \code{edina} model object.
}
