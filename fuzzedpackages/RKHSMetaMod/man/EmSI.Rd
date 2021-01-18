\name{SI_emp}
\alias{SI_emp}
\title{
Function to calculate the empirical sensitivity indices for an input or a group of inputs.
}
\description{
Calculates the empirical sensitivity indices.
}
\usage{
SI_emp(res,ErrPred)
}
\arguments{
  \item{res}{List, includes a sequence of estimated meta models, the solutions of the RKHS Ridge Group Sparse or RKHS Group Lasso problems. It should has the same format as the output of one of the functions: \code{\link[RKHSMetaMod]{pen_MetMod}}, \code{\link[RKHSMetaMod]{RKHSMetMod}} or \code{\link[RKHSMetaMod]{RKHSMetMod_qmax}}.}
  \item{ErrPred}{Matrix or NULL. If matrix, each element of the matrix is the obtained prediction error associated with one RKHS meta model in "res". It should have the same format as the output of the function \code{\link[RKHSMetaMod]{PredErr}}. Set as "NULL" by default.}
}
\details{
Details.
}
\value{
If input ErrPred\eqn{\neq}"NULL", Vector of the empirical sensitivity incdices for the meta model with the minimum Prediction error is returned. If ErrPred\eqn{=}"NULL", a list of the vectors is returned. Each vector is the obtained sensitivity indices associated with one meta model in "res".
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
\code{\link[RKHSMetaMod]{PredErr}}, \code{\link[RKHSMetaMod]{pen_MetMod}}, \code{\link[RKHSMetaMod]{RKHSMetMod}}, \code{\link[RKHSMetaMod]{RKHSMetMod_qmax}}
}
\examples{
d <- 3
n <- 50;nT <- 50
library(lhs)
X <- maximinLHS(n, d);XT <- maximinLHS(nT, d)
c <- c(0.2,0.6,0.8)
F <- 1;for (a in 1:d) F <- F*(abs(4*X[,a]-2)+c[a])/(1+c[a])
FT <- 1;for (a in 1:d) FT <- FT*(abs(4*XT[,a]-2)+c[a])/(1+c[a])
sigma <- 0.2
epsilon <- rnorm(n,0,1);Y <- F + sigma*epsilon
epsilonT <- rnorm(nT,0,1);YT <- FT + sigma*epsilonT
Dmax <- 3
kernel <- "matern"
frc <- c(10)
gamma=c(.5,.01,.001)
res <- RKHSMetMod(Y,X,kernel,Dmax,gamma,frc,FALSE)
mu <- vector()
l <- length(gamma)
for(i in 1:length(frc)){mu[i]=res[[(i-1)*l+1]]$mu}
error <- PredErr(X,XT, YT,mu,gamma, res, kernel,Dmax)
SI.minErr  <- SI_emp(res, error)
SI.minErr
SI <- SI_emp(res, NULL)
SI
}
\keyword{ Sensitivity indices }
\keyword{ Empirical }
