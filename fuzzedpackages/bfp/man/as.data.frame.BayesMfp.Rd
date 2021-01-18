\name{as.data.frame.BayesMfp}
\alias{as.data.frame.BayesMfp}

\title{Convert a BayesMfp object to a data frame}
\description{
  Convert the \code{\link{BayesMfp}} object to a data frame with the
  saved models.
}

\usage{
\method{as.data.frame}{BayesMfp}(x, row.names = NULL, ..., freq = TRUE)
}

\arguments{
  \item{x}{valid \code{\link{BayesMfp}} object}
  \item{row.names}{optional rownames (default is to keep the names of
  the \code{\link{BayesMfp}} list)}
  \item{freq}{should empirical frequencies of the models in the sampling
    path be given? (default)}
  \item{\dots}{unused}
}

\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{summary.BayesMfp}}}

\examples{
## generate a BayesMfp object
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

## get the models data frame
as.data.frame(test)
}

\keyword{methods}

