\name{plot.sim_dat}
\alias{plot.sim_dat}
\title{
Plot the true covariate effects
}
\description{
This function plots the functional form of covariate effects in four simulation scenarios.
}
\usage{
\method{plot}{sim_dat}(x, which.predictor = NULL, n.plot = 4, \dots)
}
\arguments{
  \item{x}{
an object of class "sim_dat"
}
  \item{which.predictor}{
a vector of predictor index that indicates which predictor function to plot. The vector should have integer values from 1 to p where p is the number of predictors.
}
  \item{n.plot}{
the number of predictors to be plotted (default is 4). If \code{which.predictor} is supplied, this argument is ignored.
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}
\author{
Jiacheng Wu
}

\seealso{
\code{\link{sim_dat}}
}

\examples{
#generate data
set.seed(123)
dat = sim_dat(n=100, zerof=0, scenario=1)
#plot X versus the true theta
plot.sim_dat(dat)
}
