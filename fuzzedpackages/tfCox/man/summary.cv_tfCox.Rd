\name{summary.cv_tfCox}
\alias{summary.cv_tfCox}
\title{
Summarize \code{cv_tfCox} object
}
\description{
This function summarizes \code{cv_tfCox} object and identifies the tuning parameter chosen by cross-validation.
}
\usage{
\method{summary}{cv_tfCox}(object, \dots)
}
\arguments{
  \item{object}{
an object of class "cv_tfCox"
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}

\author{
Jiacheng Wu
}

\seealso{
\code{\link{cv_tfCox}}, \code{\link{plot.cv_tfCox}}
}

\examples{
#generate data
set.seed(1234)
dat = sim_dat(n=100, zerof=0, scenario=1)

#cross-validation to choose the tuning parameter lambda with fixed alpha=1
cv = cv_tfCox(dat, ord=0, alpha=1, n.fold=2)
#summarize the cross-validation
summary(cv)
#plot the cross-validation curve
plot(cv)
}
