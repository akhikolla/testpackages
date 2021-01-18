\name{table_tennis}
\alias{table_tennis}
\alias{table_tennis_serve}
\alias{ping_pong}
\alias{profile}
\alias{profile_likelihood}
\docType{data}
\title{Match outcomes from repeated table tennis matches}
\description{Match outcomes from repeated singles table tennis matches}
\usage{data(table_tennis)}
\format{
A likelihood function corresponding to the match outcomes listed below.
}
\details{

There are four players, A, B, and C, who
play singles table tennis matches with the following results:

\itemize{
  \item A vs B, A serves, 5-1
  \item A vs B, B serves, 1-3
  \item A vs C, A serves, 4-1
  \item A vs C, C serves, 1-2
}

As discussed in vignette \code{table_tennis_serve}, we wish to assess
the importance of the serve.  The vignette presents a number of analyses
including a profile likelihood plot.

See vignette \code{table_tennis_serve} for an account of how to create
\code{table_tennis}.
  
}

\examples{
data(table_tennis)
dotchart(maxp(table_tennis))
}
\keyword{datasets}
