% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbt.ksval.R
\name{gbt.ksval}
\alias{gbt.ksval}
\title{Kolmogorov-Smirnov validation of model}
\usage{
gbt.ksval(object, y, x)
}
\arguments{
\item{object}{Object or pointer to object of class \code{ENSEMBLE}}

\item{y}{Observations to be tested}

\item{x}{design matrix for training. Must be of type \code{matrix}.}
}
\value{
Kolmogorov-Smirnov test of model
}
\description{
\code{gbt.ksval} transforms observations to U(0,1) if the model
is correct and performs a Kolmogorov-Smirnov test for uniformity.
}
\details{
Model validation of model passed as \code{object} using observations \code{y}.
Assuming the loss is a negative log-likelihood and thus a probabilistic model, 
the transformation 
\deqn{u = F_Y(y;x,\theta) \sim U(0,1),}
is usually valid. 
One parameter, \eqn{\mu=g^{-1}(f(x))}, is given by the model. Remaining parameters 
are estimated globally over feature space, assuming they are constant.
This then allow the above transformation to be exploited, so that the 
Kolmogorov-Smirnov test for uniformity can be performed.

If the response is a count model (\code{poisson} or \code{negbinom}), the transformation
\deqn{u_i = F_Y(y_i-1;x,\theta) + Uf_Y(y_i,x,\theta), ~ U \sim U(0,1)}
is used to obtain a continuous transformation to the unit interval, which, if the model is
correct, will give standard uniform random variables.
}
\examples{
## Gaussian regression:
x_tr <- as.matrix(runif(500, 0, 4))
y_tr <- rnorm(500, x_tr, 1)
x_te <- as.matrix(runif(500, 0, 4))
y_te <- rnorm(500, x_te, 1)
mod <- gbt.train(y_tr, x_tr)
gbt.ksval(mod, y_te, x_te)

}
