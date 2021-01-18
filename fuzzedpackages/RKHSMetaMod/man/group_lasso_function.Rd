\name{RKHSgrplasso}
\alias{RKHSgrplasso}
\title{
Function to fit a solution of an RKHS Group Lasso problem.
}
\description{
Fits the solution of an RKHS group lasso problem for the Gaussian regression model.
}
\usage{
RKHSgrplasso(Y, Kv, mu, maxIter, verbose)
}

\arguments{
  \item{Y}{Vector of response observations of size \eqn{n}.}
  \item{Kv}{List, includes the eigenvalues and eigenvectors of the positive definite Gram matrices \eqn{K_v, v=1,...,}vMax and their associated group names. It should have the same format as the output of the function \code{\link[RKHSMetaMod]{calc_Kv}} (see details).}
  \item{mu}{Positive scalar, value of the penalty parameter \eqn{\mu_g} in the RKHS Group Lasso problem.}
  \item{maxIter}{Integer, shows the maximum number of loops through all groups. Set as \eqn{1000} by default.}
  \item{verbose}{Logical, if TRUE, prints: the number of current iteration, active groups and convergence criterias. Set as FALSE by default.}
}
\details{
Input Kv should contain the eigenvalues and eigenvectors of positive definite Gram matrices \eqn{K_v}. It is necessary to set input correction in the function \code{\link[RKHSMetaMod]{calc_Kv}} equal to "TRUE".
}
\value{
Estimated RKHS meta model, list with \eqn{13} components:
\item{intercept}{Scalar, estimated value of intercept.}
\item{teta}{Matrix with vMax rows and \eqn{n} columns. Each row of the matrix is the estimated vector \eqn{\theta_{v}} for \eqn{v=1,...,}vMax.}
\item{fit.v}{Matrix with \eqn{n} rows and vMax columns. Each row of the matrix is the estimated value of \eqn{f_{v}=K_{v}\theta_{v}}.}
\item{fitted}{Vector of size \eqn{n}, indicates the estimator of \eqn{m}.}
\item{Norm.H}{Vector of size vMax, estimated values of the penalty norm.}
\item{supp}{Vector of active groups.}
\item{Nsupp}{Vector of the names of the active groups.}
\item{SCR}{Scalar, equals to \eqn{\Vert Y-f_{0}-\sum_{v}K_{v}\theta_{v}\Vert ^{2}}.}
\item{crit}{Scalar, indicates the value of penalized criteria.}
\item{MaxIter}{Integer, number of iterations until convergence is reached.}
\item{convergence}{TRUE or FALSE. Indicates whether the algorithm has converged or not.}
\item{RelDiffCrit}{Scalar, value of the first convergence criteria at the last iteration, \eqn{\frac{crit_{lastIter}-crit_{lastIter-1}}{crit_{lastIter-1}}}.}
\item{RelDiffPar}{Scalar, value of the second convergence criteria at the last iteration, \eqn{\Vert\frac{\theta_{lastIter}-\theta_{lastIter-1}}{\theta_{lastIter-1}}\Vert ^{2}}.}
}
\references{
Kamari, H., Huet, S. and Taupin, M.-L. (2019) RKHSMetaMod : An R package to estimate the Hoeffding decomposition of an unknown function by solving RKHS Ridge Group Sparse optimization problem. <arXiv:1905.13695>

Meier, L. Van de Geer, S. and Buhlmann, P. (2008) The group LASSO for logistic regression. Journal of the Royal Statistical Society Series B. 70. 53-71. 10.1111/j.1467-9868.2007.00627.x.
}
\author{
Halaleh Kamari
}
\note{
Note.
}

\seealso{
\code{\link[RKHSMetaMod]{calc_Kv}}
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
Kv <- calc_Kv(X, kernel, Dmax, TRUE, TRUE)
matZ <- Kv$kv
mumax <- mu_max(Y, matZ)
mug <- mumax/10
gr <- RKHSgrplasso(Y,Kv, mug , 1000, FALSE)
gr$Nsupp
}
\keyword{ Meta Model }
\keyword{RKHS Group Lasso }
\keyword{ Optimization problem }