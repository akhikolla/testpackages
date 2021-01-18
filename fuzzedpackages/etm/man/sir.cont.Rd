\name{sir.cont}
\docType{data}
\alias{sir.cont}
\title{Ventilation status in intensive care unit patients}
\description{
  Time-dependent ventilation status for intensive care unit (ICU)
  patients, a random sample from the SIR-3 study.
}
\usage{
data(sir.cont)
}
\format{
  A data frame with 1141 rows and 6 columns:
  \describe{
  \item{id:}{Randomly generated patient id}
  \item{from:}{State from which a transition occurs}
  \item{to:}{State to which a transition occurs}
  \item{time:}{Time when a transition occurs}
  \item{age:}{Age at inclusion}
  \item{sex:}{Sex. \code{F} for female and \code{M} for male}
  }
  
  The possible states are:
  
  0: No ventilation
  
  1: Ventilation
  
  2: End of stay

  And \code{cens} stands for censored observations.
}

\details{
  This data frame consists in a random sample of the SIR-3 cohort
  data. It focuses on the effect of ventilation on the length of stay
  (combined endpoint discharge/death). Ventilation status is considered
  as a transcient state in an illness-death model. 

  The data frame is directly formated to be used with the \code{etm}
  function, i.e. it is transition-oriented with one row per transition.
}

\references{
  Beyersmann, J., Gastmeier, P., Grundmann, H.,
  Baerwolff, S., Geffers, C., Behnke, M.,
  Rueden, H., and Schumacher, M. Use of multistate
  models to assess prolongation of intensive care unit stay due to
  nosocomial infection. \emph{Infection Control and Hospital
  Epidemiology}, 27:493-499, 2006.
}

\examples{
data(sir.cont)
}
\keyword{datasets}