\name{abortion}
\alias{abortion}
\docType{data}
\title{Pregnancies exposed to coumarin derivatives}
\description{
  Outcomes of pregnancies exposed to coumarin derivatives. The aim is to
  investigate whether exposition to coumarin derivatives increases the
  probability of spontaneous abortions. Apart from spontaneous abortion,
  pregnancy may end in induced abortion or live birth, leading to a
  competing risks situation. Moreover, data are left-truncated as women
  usually enter the study several weeks after conception.
}
\usage{data(abortion)}
\format{
  A data frame with 1186 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{Identification number}
    \item{\code{entry}}{Entry times into the cohort}
    \item{\code{exit}}{Event times}
    \item{\code{group}}{Group. 0: control, 1: exposed to coumarin
      derivatives}
    \item{\code{cause}}{Cause of failure. 1: induced abortion, 2: life
      birth, 3: spontaneous abortion}
  }
}
\source{
  Meiester, R. and Schaefer, C (2008). Statistical methods for
  estimating the probability of spontaneous abortion in observational
  studies -- Analyzing pregnancies exposed to coumarin
  derivatives. Reproductive Toxicology, 26, 31--35
}
\examples{
data(abortion)
}
\keyword{datasets}
