\name{summary.etmCIF}
\alias{summary.etmCIF}
\alias{print.summary.etmCIF}
\title{
  Summary function for cifETM
}
\description{
  Summary function for objects of class \code{cifETM}
}
\usage{
\S3method{summary}{etmCIF}(object, ci.fun = "cloglog",
        level = 0.95, ...)
\S3method{print}{summary.etmCIF}(x, ...)
}
\arguments{
  \item{object}{An object of class \code{etmCIF}}
  \item{ci.fun}{Transformation applied to the pointwise confidence
    intervals. On of \code{"linear", "log", "log-log",
      "cloglog"}. Default is \code{"cloglog"}.}
  \item{level}{Level of the confidence intervals. Default is 0.95.}
  \item{x}{An object of class \code{cifETM}.}
  \item{\dots}{Further arguments}
}
\value{
  A data.frame per covariate level and competing event
    \item{P}{Transition probability estimates}
  \item{var}{Variance estimates}
  \item{lower}{Lower confidence limit}
  \item{upper}{Upper confidence limit}
  \item{time}{Transition times}
  \item{n.risk}{Number of individuals at risk of experiencing a transition
    just before time \eqn{t}{t}}
  \item{n.event}{Number of events at time \eqn{t}{t}}
}
\author{
  Arthur Allignol \email{arthur.allignol@gmail.com}
}
\seealso{
  \code{\link{etmCIF}}
}
\keyword{method}
\keyword{print}
\keyword{survival}
