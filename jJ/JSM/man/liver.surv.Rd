
\name{liver.surv}
\alias{liver.surv}
\docType{data}

\title{
Prednisone versus Placebo in Liver Cirrhosis Patients - Survival Data
}

\description{
A randomized control trial in which both survival and longitudinal data were collected to examine the development of prothrombin index over time and its relationship with the survival outcome. 488 patients were randomly allocated to prednisone (251) or placebo (237) and followed until death or end of the study. \code{liver.surv} only contains the survival data of the trial, with one row per patient.
}

\format{
  A data frame with 488 observations on the following 4 variables.
  \describe{
    \item{\code{ID}}{patient ID, there are 488 patients in total.}
    \item{\code{Time}}{survival time, i.e. time to death or censoring.}
    \item{\code{death}}{death indicator: 0 denotes censoring; 1 denotes death.}
    \item{\code{Trt}}{treatment indicator with two levels: \code{placebo} and \code{prednisone}.}
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
\code{\link{liver.long}}
}

\examples{
head(liver.surv)
}

\keyword{datasets}
