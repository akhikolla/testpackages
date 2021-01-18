\name{summary.predgc}
\alias{summary.predgc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for Extracting Information from Fitted Object of Class \code{predgc}}


\description{
Output a summary data frame.
}
\usage{
\method{summary}{predgc}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{predgc} inherited from function \code{\link{predgc}}.
}

\item{...}{further arguments.}
}



\value{
A table including the following information:

  \item{pred.locs}{ prediction locations.}

  \item{predMean}{the expectation of the conditional predictive distribution.}

  \item{predCount}{ predicted counts; the closest integer that \code{predMean} rounded to.}

  \item{predVar}{ estimated variance of the prediction at prediction locations.}

 \item{predInterval.EqualTail}{ equal-tail prediction interval;
                      computed only if \code{ConfidenceLevel = TRUE}. }

 \item{predInterval.Shortest}{ shortest length prediction interval;
                      computed only if \code{ConfidenceLevel = TRUE}.}
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link{mlegc}},
\code{\link{predgc}}
}

\keyword{Prediction}
