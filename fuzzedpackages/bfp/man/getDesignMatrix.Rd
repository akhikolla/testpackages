\name{getDesignMatrix}
\alias{getDesignMatrix}

\title{Extract the design matrix of a multiple FP model.}
\description{
  Extract the (centered) design matrix of the first element (== model) of a
  \code{\link{BayesMfp}} object. This is an internal function not
  intended to be used publicly.
}
\usage{
getDesignMatrix(x, fixedColumns=TRUE, center=TRUE)
}
\arguments{
  \item{x}{a valid BayesMfp-Object of length 1 (otherwise only first
    element is recognized)}
  \item{fixedColumns}{return the fixed columns (intercept) inside the
    matrix (default) or not?}
  \item{center}{shall the non-intercept columns be centered? (default)}
}
\value{
  The design matrix with an attribute \code{shifts} containing the used
  shifts for the optional centering of the non-intercept columns. (If
  \code{center} is \code{FALSE}, the shifts vector will contain only
  zeroes.) 
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{getFpTransforms}}}
\examples{
## construct a BayesMfp object
set.seed(19)

x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)

y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 200, method="exhaustive")

## get the design matrix of the fifth best model
a <- bfp:::getDesignMatrix(test[5])
attr(a, "shifts")

## and once again but without centering
b <- bfp:::getDesignMatrix(test[5], center=FALSE)
stopifnot(all(attr(b, "shifts") == 0))
}
\keyword{utilities}
\keyword{internal}
