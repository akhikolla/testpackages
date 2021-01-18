\name{denet}
\alias{denet}
\alias{renet}
\title{
The elastic net prior distribution
}
\description{
Density function and random number generator for the elastic net prior distribution.
}
\usage{
denet(x, lambda1=1, lambda2=1, log=FALSE)

renet(n, lambda1=1, lambda2=1)
}

\arguments{
\item{x}{
vector of quantiles
}
\item{n}{
number of samples
}
\item{lambda1}{
lambda1 parameter value
}
\item{lambda2}{
lambda2 parameter value
}
\item{log}{
should the logarithm of the density be returned
}
}

\details{
The elastic net prior density has density:
\deqn{f(x)=g(\lambda_1,\lambda_2) e^[-0.5*(\lambda_1 |x| + \lambda_2 x^2)]}
}

\value{
\code{denet} gives the density of the input \code{x}. \code{renet} gives a vector of length \code{n} of random values.
}

\references{
Münch, M.M., Peeters, C.F.W., van der Vaart, A.W., and van de Wiel, M.A. (2018). Adaptive group-regularized logistic elastic net regression. arXiv:1805.00389v1 [stat.ME].
}

\author{
Magnus M. Münch <m.munch@vumc.nl>
}

\examples{
## Create data
n <- 100
x <- renet(n)
hist(x)

## Calculate density
dens <- denet(x)
plot(sort(x), dens[order(x)])
}