
\name{liver.long}
\alias{liver.long}
\docType{data}

\title{
Prednisone versus Placebo in Liver Cirrhosis Patients - Longitudinal Data
}

\description{
A randomized control trial in which both survival and longitudinal data were collected to examine the development of prothrombin index over time and its relationship with the survival outcome. 488 patients were randomly allocated to prednisone (251) or placebo (237) and followed until death or end of the study. \code{liver.long} only contains the longitudinal data of the trial, with one row per prothrombin index measurement.
}

\format{
  A data frame with 2968 observations on the following 3 variables.
  \describe{
    \item{\code{ID}}{patient ID, there are 488 patients in total.}
    \item{\code{obstime}}{time points at which the longitudinal measurements, i.e. prothrombin index, are recorded.}
    \item{\code{proth}}{prothrombin index measured at \code{obstime}.}
  }
}

\source{
Andersen, P. K., Borgan O., Gill, R. D. and Kieding, N. (1993) \emph{Statistical Models Based on Counting Processes}. New York: Springer.
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2002) Identification and efficacy of longitudinal markers for survival. \emph{Biostatistics} \bold{3}, 33--50
}

\seealso{
\code{\link{liver}},
\code{\link{liver.surv}}
}

\examples{
head(liver.long)
}

\keyword{datasets}
