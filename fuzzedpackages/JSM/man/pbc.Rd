
\name{pbc}
\alias{pbc}
\docType{data}

\title{
Mayo Clinic Primary Biliary Cirrhosis Data
}

\description{
A randomized control trial from Mayo Clinic in which both survival and longitudinal data were collected from 1974 to 1984 to study the progression of primary biliary cirrhosis.
}

\format{
  A data frame with 1945 observations on the following 16 variables.
  \describe{
    \item{\code{ID}}{patient ID, there are 312 patients in total.}
    \item{\code{Time}}{survival time (in years), i.e. time to death, transplantion or censoring.}
    \item{\code{death}}{death indicator: 0 denotes transplantion or censoring; 1 denotes death.}
    \item{\code{obstime}}{time points at which the longitudinal measurements, e.g. serum bilirubin, albumin and alkaline phosphatase, are recorded.}
    \item{\code{serBilir}}{serum bilirubin measured at \code{obstime} (mg/dl).}
    \item{\code{albumin}}{albumin measured at \code{obstime} (gm/dl).}
    \item{\code{alkaline}}{alkaline phosphatase measured at \code{obstime} (U/litter).}
    \item{\code{platelets}}{platelets per cubic measured at \code{obstime} (ml/1000).}
    \item{\code{drug}}{drug indicator with two levels: \code{placebo} and \code{D-penicil}.}
    \item{\code{age}}{age of patient at study entry.}
    \item{\code{gender}}{gender indicator with two levels: \code{male} and \code{female}.}
    \item{\code{ascites}}{ascites indicator with two levels: \code{No} and \code{Yes}.}
    \item{\code{hepatom}}{hepatomegaly indicator with two levels: \code{No} and \code{Yes}.}
    \item{\code{start}}{same with \code{obstime}, starting time of the interval which contains the time of the logitudinal measurements.}
    \item{\code{stop}}{ending time of the interval which contains the time of the longitudinal measurements.}
    \item{\code{event}}{event indicator suggesting whether the event-of-interest, i.e. death, happens in the interval given by \code{start} and \code{stop}.}
  }
}

\source{
\url{http://lib.stat.cmu.edu/datasets/pbcseq}

Fleming, T. and Harrington, D. (1991) \emph{Counting Processes and Survival Analysis}. Wiley, New York.
}

\references{
Murtaugh, P. A., Dickson, E. R., Van Dam, G. M., Malincho, M., Grambsch, P. M., Langworthy, A. L., and Gips, C. H. (1994) Primary biliary cirrhosis: Prediction of short-term survival based on repeated patient visits. \emph{Hepatology} \bold{20}, 126--134.

Therneau, T. and Grambsch, P. (2000) \emph{Modeling Survival Data: Extending the Cox Model}. New York: Springer.

Ding, J. and Wang, J. L. (2008) Modeling longitudinal data with nonparametric multiplicative random effects jointly with survival data. \emph{Biometrics} \bold{64}, 546--556.
}

\examples{
head(pbc)
}

\keyword{datasets}
