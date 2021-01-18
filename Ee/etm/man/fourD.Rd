\name{fourD}
\alias{fourD}
\docType{data}
\title{
  Placebo data from the 4D study
}
\description{
  Data from the placebo group of the 4D study. This study aimed at
  comparing atorvastatin to placebo for patients with type 2 diabetes
  and receiving hemodialysis in terms of cariovascular events. The
  primary endpoint was a composite of
  death from cardiac causes, stroke and non-fatal myocardial infarction.
  Competing event was death from other causes.
}
\usage{data(fourD)}
\format{
  A data frame with 636 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{Patients' id number}
    \item{\code{sex}}{Patients' gender}
    \item{\code{age}}{Patients' age}
    \item{\code{medication}}{Character vector indicating treatment
      affiliation. Here only equal to \code{"Placebo"}}
    \item{\code{status}}{Status at the end of the follow-up. 1 for the
      event of interest, 2 for death from other causes and 0 for
      censored observations}
    \item{\code{time}}{Survival time}
    \item{\code{treated}}{Numeric vector indicated whether patients are
      treated or not. Here always equal to zero}
  }
}

\source{
  Wanner, C., Krane, V., Maerz, W., Olschewski, M., Mann, J., Ruf, G.,
  Ritz, E (2005). Atorvastatin in patients with type 2 diabetes mellitus
  undergoing hemodialysis. New England Journal of Medicine, 353(3),
  238--248.
}
\references{
  Allignol, A., Schumacher, M., Wanner, C., Dreschler, C. and
  Beyersmann, J. (2010). Understanding competing risks: a simulation
  point of view. Research report.
}
\examples{
data(fourD)
}
\keyword{datasets}
