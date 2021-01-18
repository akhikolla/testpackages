\name{plot.cv_tfCox}
\alias{plot.cv_tfCox}
\title{
Plots Cross-Validation Curve for Object of Class "cv_tfCox"
}
\description{
This function plots the cross-validation curve for models fitted by a range of tuning parameter lambda using \code{\link{cv_tfCox}}. The cross-validation error with +/-1 standard error is plotted for each value of lambda. The dotted vertical line indicates the chosen lambda corresponding to the minimum cross-validation error.
}
\usage{
\method{plot}{cv_tfCox}(x, showSE=F, \dots)
}
\arguments{
  \item{x}{
an object of class "cv_tfCox"
}
  \item{showSE}{
a logical (TRUE or FALSE) for whether the standard errors of the curve should be plotted
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}


\author{
Jiacheng Wu
}

\seealso{
\code{\link{cv_tfCox}}
}
\examples{
#generate data
set.seed(123)
dat = sim_dat(n=100, zerof=0, scenario=1)

#fit piecewise constant functions
#cross-validation to choose the tuning parameter lambda with fixed alpha=1
cv = cv_tfCox(dat, ord=0, alpha=1, n.fold=2, seed=123)
plot(cv, showSE=TRUE)
}
