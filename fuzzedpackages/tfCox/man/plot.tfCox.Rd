\name{plot.tfCox}
\alias{plot.tfCox}
\title{
Plot Fitted Functions from Class "tfCox"
}
\description{
This function plots the fitted functions from a model estimated by \code{\link{tfCox}}.
}
\usage{
\method{plot}{tfCox}(x, which.lambda=1, which.predictor = NULL, n.plot = 4, \dots)
}
\arguments{
  \item{x}{
an object of class "tfCox"
}
  \item{which.lambda}{
the index for the model of interest to be plotted. \code{which.lambda} corresponds to the model fit in \code{lambda.seq} and should be integer between 1 to \code{length(fit$lambda.seq)}. In other words, the fit from \code{fit$theta.list[[which.lambda]]} will be plotted.
}
  \item{which.predictor}{
a vector of predictor index that indicates which predictor function to plot. The vector should have integer values from 1 to p where p is the number of predictors.
}
  \item{n.plot}{
the number of predictors to be plotted (default is 4). Note that only those non-zero estimated functions will be plotted. If \code{which.predictor} is supplied, this argument is ignored.
  }
    \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}

\author{
Jiacheng Wu
}

\seealso{
\code{\link{tfCox}}
}
\examples{
#generate data
set.seed(123)
dat = sim_dat(n=100, zerof=0, scenario=1)

fit = tfCox(dat, ord=0, alpha=1, lambda.seq=0.04)
plot(fit, n.plot=4)

}
