
\name{liver}
\alias{liver}
\docType{data}

\title{
Prednisone versus Placebo in Liver Cirrhosis Patients
}

\description{
A randomized control trial in which both survival and longitudinal data were collected to examine the development of prothrombin index over time and its relationship with the survival outcome. 488 patients were randomly allocated to prednisone (251) or placebo (237) and followed until death or end of the study.
}

\format{
  A data frame with 2968 observations on the following 9 variables.
  \describe{
    \item{\code{ID}}{patient ID, there are 488 patients in total.}
    \item{\code{Time}}{survival time, i.e. time to death or censoring.}
    \item{\code{death}}{death indicator: 0 denotes censoring; 1 denotes death.}
    \item{\code{obstime}}{time points at which the longitudinal measurements, i.e. prothrombin index, are recorded.}
    \item{\code{proth}}{prothrombin index measured at \code{obstime}.}
    \item{\code{Trt}}{treatment indicator with two levels: \code{placebo} and \code{prednisone}.}
    \item{\code{start}}{same with \code{obstime}, starting time of the interval which contains the time of the prothrombin index measurement.}
    \item{\code{stop}}{ending time of the interval which contains the time of the prothrombin index measurement.}
    \item{\code{event}}{event indicator suggesting whether the event-of-interest, i.e. death, happens in the interval given by \code{start} and \code{stop}.}
  }
}

\source{
Andersen, P. K., Borgan O., Gill, R. D. and Kieding, N. (1993) \emph{Statistical Models Based on Counting Processes}. New York: Springer.
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2002) Identification and efficacy of longitudinal markers for survival. \emph{Biostatistics} \bold{3}, 33--50
}

\seealso{
\code{\link{liver.long}},
\code{\link{liver.surv}}
}

\examples{
head(liver)
}

\keyword{datasets}
