\name{pen_MetMod}
\alias{pen_MetMod}

\title{
Function to fit a solution of the RKHS Ridge Group Sparse problem.
}
\description{
Fits the solution of the RKHS ridge group sparse optimization problem for the Gaussian regression model.}
\usage{
pen_MetMod(Y, Kv, gamma, mu, resg, gama_v, mu_v, maxIter, verbose, calcStwo)
}
\arguments{
  \item{Y}{Vector of response observations of size \eqn{n}.}
  \item{Kv}{List, includes the eigenvalues and eigenvectors of the positive definite Gram matrices \eqn{K_v, v=1,...,}vMax and their associated group names. It should have the same format as the output of the function \code{\link[RKHSMetaMod]{calc_Kv}} (see details).}
  \item{gamma}{Vector of positive scalars. Values of the penalty parameter \eqn{\gamma} in decreasing order.}
  \item{mu}{Vector of positive scalars. Values of the penalty parameter \eqn{\mu} in decreasing order.}
  \item{resg}{List of initial parameters, includes the \code{\link[RKHSMetaMod]{RKHSgrplasso}} objects for each value of the penalty parameter \eqn{\mu}.}
  \item{gama_v}{Scalar zero or vector of vMax positive scalars, considered as weights for the Ridge penalty. Set to zero, to consider no weights, i.e. all weights equal to \eqn{1}.}
  \item{mu_v}{Scalar zero or a vector with vMax scalars, considered as weigths of Sparse Group penalty. Set to zero, to consider no weights, i.e. all weights equal to \eqn{1}.}
  \item{maxIter}{Integer, shows the maximum number of loops through initial active groups at the first step and maximum number of loops through all groups at the second step. Set as \eqn{1000} by default.}
  \item{verbose}{Logical, if TRUE, for each pair of penalty parameters \eqn{(\mu,\gamma)} it prints: the number of current iteration, active groups and convergence criterias. Set as FALSE by default.}
  \item{calcStwo}{Logical, if TRUE, the program does a second step after convergence: the algorithm is done over all groups by taking the estimated parameters at the first step as initial values. Set as FALSE by default.}
}
\details{
Input Kv should contain the eigenvalues and eigenvectors of positive definite Gram matrices \eqn{K_v}. It is necessary to set input "correction" in the function \code{\link[RKHSMetaMod]{calc_Kv}} equal to "TRUE".
}
\value{
List of l components, with l equals to the number of pairs of the penalty parameters \eqn{(\mu,\gamma)}. Each component of the list is a list of \eqn{3} components "mu", "gamma" and "Meta-Model":
\item{mu}{Positive scalar, an element of the input vector mu associated with the estimated Meta-Model.}
\item{gamma}{Positive scalar, an element of the input vector gamma associated with the estimated Meta-Model.}
\item{Meta-Model}{Estimated meta model associated with penalty parameters mu and gamma. List of \eqn{16} components:}
\item{intercept}{Scalar, estimated value of intercept.}
\item{teta}{Matrix with vMax rows and \eqn{n} columns. Each row of the matrix is the estimated vector \eqn{\theta_{v}} for \eqn{v=1,...,}vMax.}
\item{fit.v}{Matrix with \eqn{n} rows and vMax columns. Each row of the matrix is the estimated value of \eqn{f_{v}=K_{v}\theta_{v}}.}
\item{fitted}{Vector of size \eqn{n}, indicates the estimator of \eqn{m}.}
\item{Norm.n}{Vector of size vMax, estimated values for the Ridge penalty norm.}
\item{Norm.H}{Vector of size vMax, estimated values for the Group Sparse penalty norm.}
\item{supp}{Vector of active groups.}
\item{Nsupp}{Vector of the names of the active groups.}
\item{SCR}{Scalar equals to \eqn{\Vert Y-f_{0}-\sum_{v}K_{v}\theta_{v}\Vert ^{2}}.}
\item{crit}{Scalar indicates the value of the penalized criteria.}
\item{gamma.v}{Vector of size vMax, coefficients of the Ridge penalty norm, \eqn{\sqrt{n}\gamma\times}gama_v.}
\item{mu.v}{Vector of size vMax, coefficients of the Group Sparse penalty norm, \eqn{n\mu\times}mu_v.}
\item{iter}{List of three components if calcStwo\eqn{=}TRUE \eqn{(}two components if calcStwo\eqn{=}FALSE\eqn{)}: maxIter, number of iterations until convergence is reached at first step and the number of iterations until convergence is reached at second step \eqn{(}maxIter, and the number of iterations until convergence is reached at first step\eqn{)}.}
\item{convergence}{TRUE or FALSE. Indicates whether the algorithm has converged or not.}
\item{RelDiffCrit}{List of two components if calcStwo\eqn{=}TRUE \eqn{(}one component if calcStwo\eqn{=}FALSE\eqn{)}: value of convergence criteria at the last iteration of each step, \eqn{\Vert\frac{\theta_{lastIter}-\theta_{lastIter-1}}{\theta_{lastIter-1}}\Vert ^{2}}.}
\item{RelDiffPar}{List of two components if calcStwo\eqn{=}TRUE (one component if calcStwo\eqn{=}FALSE): value of convergence criteria at the last iteration, \eqn{\frac{crit_{lastIter}-crit_{lastIter-1}}{crit_{lastIter-1}}} of each step.}
}
\references{
Huet, S. and Taupin, M. L. (2017) Metamodel construction for sensitivity analysis. ESAIM:
Procs 60, 27-69.

Kamari, H., Huet, S. and Taupin, M.-L. (2019) RKHSMetaMod : An R package to estimate the Hoeffding decomposition of an unknown function by solving RKHS Ridge Group Sparse optimization problem. <arXiv:1905.13695>
}
\author{
Halaleh Kamari
}
\note{
Note.
}

\seealso{
\code{\link[RKHSMetaMod]{calc_Kv}}, \code{\link[RKHSMetaMod]{RKHSgrplasso}}
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
Kv <- calc_Kv(X, kernel, Dmax, TRUE,TRUE, tol = 1e-08)
vMax <- length(Kv$names.Grp)
matZ <- Kv$kv
mumax <- mu_max(Y, matZ)
mug1 <- mumax/10
mug2 <- mumax/100
gr1 <- RKHSgrplasso(Y,Kv, mug1)
gr2 <- RKHSgrplasso(Y,Kv, mug2)
gamma <- c(.5,.01,.001)
#rescaling the penalty parameter
mu <- c(mug1/sqrt(n),mug2/sqrt(n))
resg <- list(gr1,gr2)
res <- pen_MetMod(Y,Kv,gamma,mu,resg,0,0)
l <- length(res)
for(i in 1:l){print(res[[i]]$mu)}
for(i in 1:l){print(res[[i]]$gamma)}
for(i in 1:l){print(res[[i]]$`Meta-Model`$Nsupp)}
gama_v <- rep(1,vMax)
mu_v <- rep(1,vMax)
res.w <- pen_MetMod(Y,Kv,gamma,mu,resg,gama_v,mu_v)
for(i in 1:l){print(res.w[[i]]$`Meta-Model`$Nsupp)}
}
\keyword{ Meta Model }
\keyword{ RKHS Ridge Group Sparse}
\keyword{ Optimization problem }
