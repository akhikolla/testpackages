
\name{epilepsy}
\alias{epilepsy}
\docType{data}

\title{
CBZ versus LTG in Epilepsy Patients
}

\description{
A randomised control trial, the SANAD (standard and new antiepileptic drugs) study, in which both survival and longitudinal data were collected to investigate the effect of drug titration on the relative effects of two antiepileptic drugs, namely \code{CBZ} (carbamazepine, a standard drug) and \code{LTG} (lamotrigine, a new drug), on treatment failure. Treatment failure, i.e. withdrawal of the randomized drug, is the event of interest. Two main reasons for withdrawal are unacceptable adverse effects (UAE) and inadequate seizure control (ISC).
}

\format{
  A data frame with 2797 observations on the following 16 variables.
  \describe{
    \item{\code{ID}}{patient ID, there are 605 patients in total.}
    \item{\code{Time}}{survival time, i.e. time to withdrawal or censoring.}
    \item{\code{withdraw}}{withdrawal indicator: 0 denotes censoring; 1 denotes withdrawal.}
    \item{\code{withdrawUAE}}{withdrawal due to UAE indicator: 1 denotes withdrawal due to UAE; 0 otherwise.}
    \item{\code{withdrawISC}}{withdrawal due to ISC indicator: 1 denotes withdrawal due to ISC; 0 otherwise.}
    \item{\code{obstime}}{time points at which the longitudinal measurements, i.e. the dose, are recorded.}
    \item{\code{dose}}{calibrated dose measured at \code{obstime}.}
    \item{\code{drug}}{drug indicator with two levels: \code{CBZ} and \code{LTG}.}
    \item{\code{age}}{age of patient at study entry.}
    \item{\code{gender}}{gender indicator with two levels: \code{male} and \code{female}.}
    \item{\code{disab}}{learning disability indicator with two levels: \code{No} and \code{Yes}.}
    \item{\code{start}}{same with \code{obstime}, starting time of the interval which contains the time of the dose measurement.}
    \item{\code{stop}}{ending time of the interval which contains the time of the dose measurement.}
    \item{\code{event}}{event indicator suggesting whether the event-of-interest, i.e. withdrawal, happens in the interval given by \code{start} and \code{stop}.}
    \item{\code{eventUAE}}{event indicator suggesting whether the event-of-interest, i.e. withdrawal due to UAE, happens in the interval given by \code{start} and \code{stop}.}
    \item{\code{eventISC}}{event indicator suggesting whether the event-of-interest, i.e. withdrawal due to ISC, happens in the interval given by \code{start} and \code{stop}.}
    }
}

\source{
Marson, A. G., AI-Kharusi, A. M., Alwaidh, M., Appleton, R., Baker, G. A., Chadwick, D. W., Cramp, C., Cockerell, O. C. Cooper, P. N., Doughty, J. et al. (2007) The SANAD study of effectiveness of carbamazepine, gabapentin, lamotrigine, oxcarbazepine, or topiramate for treatment of partial epilepsy: an unblinded randomised controlled trial. \emph{The Lancet} \bold{369}, 1000--1015.
}

\references{
Williamson P. R., Kolamunnage-Dona R., Philipson P. and Marson A. G. (2008) Joint modelling of longitudinal and competing risks data. \emph{Statistics in Medicine} \bold{27}, 6426--6438.
}

\examples{
head(epilepsy)
}

\keyword{datasets}
