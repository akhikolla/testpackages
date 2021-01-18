\name{etmCIF}
\alias{etmCIF}
\title{
  Cumulative incidence functions of competing risks
}
\description{
  \code{etmCIF} is a wrapper around the \code{etm} function for
  facilitating the computation of the cumulative incidence functions in
  the competing risks framework.
}
\usage{
etmCIF(formula, data, etype, subset, na.action, failcode = 1)
}
\arguments{
  \item{formula}{A \code{formula} object, that must have a \code{Surv}
    object on the left of ~ operator, and a discrete covariate (or 1) on
    the right. The status indicator should be 1 (or TRUE) for an event
    (whatever the type of this event, 0 (or FALSE) for censored
    observations.)}
  \item{data}{A data.frame in which to interpret the terms of the
    formula}
  \item{etype}{Competing risks event indicator. When the status
    indicator is 1 (or TRUE) in the formula, \code{etype} describes the
    type of event, otherwise, for censored observation, the value of
    \code{etype} is ignored}
  \item{subset}{Expression saying that only a subset of the data should
    be used.}
  \item{na.action}{Missing-data filter function. Default is
    \code{options()$na.action}.}
  \item{failcode}{Indicates the failure type of interest. Default is
    one. This option is only relevant for some options of the
    \code{plot} function.}
}
\details{
    This function computes the cumulative incidence functions in a
  competing risks setting using the \code{etm} machinery, without having
  to specify the matrix of possible transitions and using the more usual
  formula specification with \code{Surv}
}
\value{
  Returns a list of \code{etm} objects (1 per covariate level) plus
  additional informations:
  \item{failcode}{As in function call}
  \item{call}{Function call}
  \item{X}{A matrix giving the name of the covariate (if present) and
    the levels of this covariate.}
}
\author{
  Arthur Allignol \email{arthur.alignol@gmail.com}
}
\seealso{
  \code{\link{etm}}, \code{\link{print.etmCIF}},
  \code{\link{summary.etmCIF}}, \code{\link{plot.etmCIF}}
}
\examples{

data(abortion)

cif.ab <- etmCIF(survival::Surv(entry, exit, cause != 0) ~ group, abortion,
                 etype = cause, failcode = 3)

cif.ab

plot(cif.ab, ci.type = "bars", pos.ci = 24,
     col = c(1, 2), lty = 1, curvlab = c("Control", "Exposed"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}