\name{trim.best}
\alias{trim.best}
\title{
Trim and refit best model chosen from grid
}
\description{
Get the best fitting model (based on minimum BIC) from grid of models,
trim off mediators with small effects (|alpha*beta|  < mediator.epsilon), then
refit with specified lambda penalty, keeping all other parameters fixed
at values specified when grid object was created.
}
\usage{
trim.best(obj, lambda = 0.0, mediator.epsilon = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
regmed.grid object created by regmed.grid
}
\item{lambda}{
  lambda penalty parameter
}
  \item{mediator.epsilon}{
tolerance to remove mediator if |alpha*beta| < mediator.epsilon
}
}
\value{
An object of regmed class
}

\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}

\seealso{
\code{\link{regmed.grid}}
}
\examples{
data(regmed_example)

y <- regmed_example$y
x <- regmed_example$x
med <- regmed_example[, -c(1,2)]
fit.grid <- regmed.grid(x, med, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)
fit.trim <- trim.best(fit.grid)
data(regmed_example)

y <- regmed_example$y
x <- regmed_example$x
med <- regmed_example[, -c(1,2)]
fit.grid <- regmed.grid(x, med, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)
fit.trim <- trim.best(fit.grid)
summary(fit.trim)
}

\keyword{ models }
