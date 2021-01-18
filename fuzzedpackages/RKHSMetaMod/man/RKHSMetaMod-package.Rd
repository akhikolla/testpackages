\name{RKHSMetaMod-package}
\alias{RKHSMetaMod-package}
\alias{RKHSMetaMod}
\docType{package}
\title{
Set of Rcpp and R functions to produce a sequence of meta models that are the solutions of the RKHS Ridge Group Sparse or RKHS Group Lasso optimization problems, calulate their associated prediction errors as well as their empirical sensitivity indices.
}
\description{
Fits a meta model to an unknown model \eqn{m} by solving the ridge group sparse (or group lasso) optimization problem based on the reproducing kernel Hilbert spaces (RKHS), for the Gaussian regression model : \deqn{Y=m(X)+\sigma\varepsilon,}
where variables \eqn{X_{1},...,X_{d}} are independent and uniformly distributed on \eqn{[0,1]} and are independent of \eqn{\varepsilon}'s.
  
We define the ridge group sparse criteria by :
\deqn{C(f_{0},\theta)=\Vert Y-f_{0}I_{n}-\sum_{v\in\mathcal{P}}K_{v}\theta_{v}\Vert ^2+\sqrt{n}\gamma\sum_{v\in\mathcal{P}}\Vert K_v\theta_{v}\Vert +n\mu\sum_{v\in\mathcal{P}}\Vert K_{v}^{1/2}\theta_{v}\Vert,}
and the group lasso criteria is obtained by setting \eqn{\gamma=0} in the criteria above. We set \eqn{\mu_g=\sqrt{n}\mu} to be the group lasso penalty parameter.

For each pair of the penalty parameters \eqn{(\mu,\gamma)} in the ridge group sparse criteria, one meta model, called RKHS meta model, is calculated. The RKHS meta model is an additive model with at most vMax groups. It satisfies the properties of the Hoeffding decomposition, and its terms estimate the terms in the Hoeffding decomposition of the function \eqn{m}.

These estimators are evaluated using a testing dataset. That is, the prediction error is calculated for each RKHS meta model and the one with the minimum prediction error is the "best" estimator for the true model \eqn{m}. It provides a function that estimates the empirical sensitivity indices of the "best" RKHS meta model as an approximation of the true sensitivity indices.
}
\details{
Details.
}
\author{
Halaleh Kamari.

Maintainer: \email{halaleh.kamari@univ-evry.fr}
}
\references{
Kamari, H., Huet, S. and Taupin, M.-L. (2019) RKHSMetaMod : An R package to estimate the Hoeffding decomposition of an unknown function by solving RKHS Ridge Group Sparse optimization problem. <arXiv:1905.13695>
}
\keyword{ Package }
\keyword{ Gaussian Regresson model }
\keyword{ Hoeffding decomposition }
\keyword{ Optimization problem }
\keyword{ Meta Model }
\keyword{ Prediction error }
\keyword{ Sensitivity indices }
\keyword{ RKHS Ridge Group Sparse}
\keyword{ RKHS Group Lasso}
\seealso{
\code{\link[RKHSMetaMod]{RKHSMetaMod}}
}
\examples{
d <- 3
n <- 50
library(lhs)
X <- maximinLHS(n, d)
c <- c(0.2,0.6,0.8)
F <- 1;for (a in 1:d) F <- F*(abs(4*X[,a]-2)+c[a])/(1+c[a])
epsilon <- rnorm(n,0,1);sigma <- 0.2
Y <- F + sigma*epsilon
Dmax <- 3
kernel <- "matern"
frc <- c(10,100)
gamma <- c(.5,.01,.001,0)
result <- RKHSMetMod(Y,X,kernel,Dmax,gamma,frc,FALSE)  
}
