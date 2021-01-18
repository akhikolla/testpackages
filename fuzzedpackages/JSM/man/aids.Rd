
\name{aids}
\alias{aids}
\docType{data}

\title{ddI versus ddC in HIV-infected Patients}

\description{
  A randomized clinical trial in which both survival and longitudinal data were collected to compare the efficacy and 
  safety of two antiretroviral drugs, namely \code{ddI} (didanosine) and \code{ddC} (zalcitabine), in treating HIV-infected patients intolerant or failing zidovudine (AZT) therapy.
}

\format{
  A data frame with 1405 observations on the following 12 variables.
  \describe{
    \item{\code{ID}}{patient ID, there are 467 patients in total.}
    \item{\code{Time}}{survival time, i.e. time to death or censoring.}
    \item{\code{death}}{death indicator: 0 denotes censoring; 1 denotes death.}
    \item{\code{obstime}}{time points at which the longitudinal measurements, i.e. CD4 cell counts, are recorded.}
    \item{\code{CD4}}{CD4 cell counts measured at \code{obstime}.}
    \item{\code{drug}}{drug indicator with two levels: \code{ddI} and \code{ddC}.}
    \item{\code{gender}}{gender indicator with two levels: \code{male} and \code{female}.}
    \item{\code{prevDiag}}{AIDS diagnosis at study entry indicator with two levels: \code{AIDS} and \code{noAIDS}.}
    \item{\code{AZT}}{AZT intolerance/failure indicator with two levels: \code{intolerance} and \code{failure}.}
    \item{\code{start}}{same with \code{obstime}, starting time of the interval which contains the time of the CD4 cell count measurement.}
    \item{\code{stop}}{ending time of the interval which contains the time of the CD4 cell count measurement.}
    \item{\code{event}}{event indicator suggesting whether the event-of-interest, i.e. death, happens in the interval given by \code{start} and \code{stop}.}
  }
}

\source{
Goldman, A., Carlin, B., Crane, L., Launer, C., Korvick, J., Deyton, L. and Abrams, D. (1996) Response of CD4+ and
clinical consequences to treatment using ddI or ddC in patients with advanced HIV infection. \emph{Journal of Acquired
Immune Deficiency Syndromes and Human Retrovirology} \bold{11}, 161--169.
}

\references{
Guo, X. and Carlin, B. (2004) Separate and joint modeling of longitudinal and event time data using standard
computer packages. \emph{The American Statistician} \bold{58}, 16--24.

Xu, C., Baines, P. D. and Wang, J. L. (2014) Standard error estimation using the EM algorithm for the joint modeling of survival and longitudinal data. \emph{Biostatistics} \bold{15}, 731--744
}

\examples{
head(aids)
}

\keyword{datasets}
