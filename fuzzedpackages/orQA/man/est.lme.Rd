\name{est.lme}
\alias{est.lme}
\title{
Estimate order restricted variance components for a two way mixed model
with interaction. 
}
\description{
Estimates variance components of a two way cross classification mixed
model with an order restricted fixed effect, a random effect and random
interaction term. 
}
\usage{
est.lme(y,ia,ib)
}
\arguments{
  \item{y}{
    a numeric matrix, for the lines of which univariate variance
    decompositions should be estimated
  }
  \item{ia}{
    integer vector of length ncol(y) specifying the ordered levels of
    the order restricted fixed effect as integers between 1 and n where
    n is the number of levels
  }
  \item{ib}{
    integer vector or factor specifying the levels of the random term
    (e.g. individuals)
  }
}
\details{
est.lme estimates the variance components for a two way cross
classification mixed model with random interaction. The order
restriction on the fixed effect used to improve the estimates by pooling
levels of the fixed effect using isotonic regression according to the
observed order in the measurements. Estimation is then done using
functionality provided by the package <nlme>.
}
\value{
  \item{sb}{
    estimates of the random effect variance component}
  \item{sg}{
    estimates of the interaction term variance component}
  \item{se}{
    estimates of the residual error variance}
  \item{tss}{
    total sum of squares (useful to obtain normalized scales between
    estimates from different measurements)}
}
\references{
  Klinglmueller, F., Tuechler, T., Posch, M. (2010)
  "Cross Platform Comparison Of Microarray Data Using Order
    Restricted Inference" Under Review

  Pinheiro, J., Bates, D., DebRoy, S., Sarkar, D.,R Development Core
    Team (2010)
    "nlme: Linear and Nonlinear Mixed Effects Models"

  Barlow, R. E., Bartholomew, D. J., Bremner, J. M., and Brunk, H. D. (1972)
    "Statistical inference under order restrictions"; Wiley, London.

  Robertson, T., Wright,F. T. and Dykstra, R. L. (1988)
    "Order Restricted Statistical Inference"; Wiley, New York.
}
\author{
  Florian Klinglmueller <float_at_lefant.net>
}
\examples{
g <- rep(1:4,each=10)
r <- rep(rep(1:2,each=5),4)

## No differences (global null)
nulldata <- matrix(rnorm(400),nc=40)

## estimation
res <- est.lme(nulldata,g,r)
round(apply(res,2,summary),2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

