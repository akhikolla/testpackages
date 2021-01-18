\name{PredErr}
\alias{PredErr}
\title{
Function to calculate the prediction error.
}
\description{
Computes the prediction error by considering a testing dataset.
}
\usage{
PredErr(X, XT, YT, mu, gamma, res, kernel, Dmax)
}
\arguments{
  \item{X}{Matrix of observations with \eqn{n} rows and \eqn{d} columns.}
  \item{XT}{Matrix of observations of the testing dataset with \eqn{n^{test}} rows and \eqn{d} columns.}
  \item{YT}{Vector of response observations of testing dataset of size \eqn{n^{test}}.}
  \item{mu}{Vector of positive scalars. Values of the Group Sparse penalty parameter in decreasing order. See function \code{\link[RKHSMetaMod]{RKHSMetMod}}.}
  \item{gamma}{Vector of positive scalars. Values of the Ridge penalty parameter in decreasing order. See function \code{\link[RKHSMetaMod]{RKHSMetMod}}.}
  \item{res}{List, includes a squence of estimated meta models for the learning dataset, using RKHS Ridge Group Sparse or RKHS Group Lasso algorithm, associated with the penalty parameters mu and gamma. It should have the same format as the output of one of the functions: \code{\link[RKHSMetaMod]{pen_MetMod}}, \code{\link[RKHSMetaMod]{RKHSMetMod}} or \code{\link[RKHSMetaMod]{RKHSMetMod_qmax}}.}
  \item{kernel}{Character, shows the type of the reproducing kernel: matern, brownian, gaussian, linear, quad. The same kernel should be chosen as the one used for the learning dataset. See function \code{\link[RKHSMetaMod]{calc_Kv}}.}
  \item{Dmax}{Integer between \eqn{1} and \eqn{d}. The same Dmax should be chosen as the one used for learning dataset. See function \code{\link[RKHSMetaMod]{calc_Kv}}.}
}
\details{
Details.
}
\value{
Matrix of the prediction errors is returned. Each element of the matrix is the obtained prediction error associated with one RKHS meta model in "res".
}
\references{
Kamari, H., Huet, S. and Taupin, M.-L. (2019) RKHSMetaMod : An R package to estimate the Hoeffding decomposition of an unknown function by solving RKHS Ridge Group Sparse optimization problem. <arXiv:1905.13695>
}
\author{
Halaleh Kamari
}
\note{
Note.
}

\seealso{
\code{\link[RKHSMetaMod]{calc_Kv}}, \code{\link[RKHSMetaMod]{pen_MetMod}}, \code{\link[RKHSMetaMod]{RKHSMetMod}}, \code{\link[RKHSMetaMod]{RKHSMetMod_qmax}}
}
\examples{
d <- 3
n <- 50
nT <- 50
library(lhs)
X <- maximinLHS(n, d)
XT <- maximinLHS(nT, d)
c <- c(0.2,0.6,0.8)
F <- 1;for (a in 1:d) F <- F*(abs(4*X[,a]-2)+c[a])/(1+c[a])
FT <- 1;for (a in 1:d) FT <- FT*(abs(4*XT[,a]-2)+c[a])/(1+c[a])
sigma <- 0.2
epsilon <- rnorm(n,0,1);Y <- F + sigma*epsilon
epsilonT <- rnorm(nT,0,1);YT <- FT + sigma*epsilonT
Dmax <- 3
kernel <- "matern"
frc <- c(10,100)
gamma <- c(.5,.01,.001)
res <- RKHSMetMod(Y,X,kernel,Dmax,gamma,frc,FALSE)
mu <- vector()
l <- length(gamma)
for(i in 1:length(frc)){mu[i]=res[[(i-1)*l+1]]$mu}
error <- PredErr(X,XT, YT,mu,gamma, res, kernel,Dmax)
error
}

\keyword{ Prediction error }
\keyword{ Testing data set }
