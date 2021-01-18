\name{grplasso_q}
\alias{grplasso_q}
\title{
Function to fit a solution with q active groups of an RKHS Group Lasso problem.
}
\description{
Fits a solution of the group lasso problem based on RKHS, with \eqn{q} active groups in the obtained solution for the Gaussian regression model. It determines \eqn{\mu_{g}(q)}, for which the number of active groups in the solution of the RKHS group lasso problem is equal to \eqn{q}, and returns the RKHS meta model associated with \eqn{\mu_{g}(q)}.
}
\usage{
grplasso_q(Y, Kv, q, rat, Num)
}

\arguments{
  \item{Y}{Vector of response observations of size \eqn{n}.}
\item{Kv}{List of eigenvalues and eigenvectors of positive definite Gram matrices \eqn{K_v} and their associated group names. It should have the same format as the output of the function \code{\link[RKHSMetaMod]{calc_Kv}} (see details).}
\item{q}{Integer, the number of active groups in the obtained solution.}
\item{rat}{Positive scalar, used to restrict the minimum value of \eqn{\mu_g}, to be evaluted in the RKHS Group Lasso algorithm, \eqn{\mu_{min}=\mu_{max}/rat}. The value \eqn{\mu_{max}} is calculated inside the program, see function \code{\link[RKHSMetaMod]{mu_max}}.}
  \item{Num}{Integer, used to restrict the number of different values of the penalty parameter \eqn{\mu_g} to be evaluated in the RKHS Group Lasso algorithm, until it achieves \eqn{\mu_g(q)}: for Num \eqn{= 1} the program is done for \eqn{3} values of \eqn{\mu_g}, \eqn{\mu_{1}=(\mu_{min}+\mu_{max})/2}, \eqn{\mu_{2}=(\mu_{min}+\mu_{1})/2} or \eqn{\mu_{2}=(\mu_{1}+\mu_{max})/2} depending on the value of \eqn{q} associated with \eqn{\mu_{1}}, \eqn{\mu_{3}=\mu_{min}}.}
}
\details{
Input Kv should contain the eigenvalues and eigenvectors of positive definite Gram matrices \eqn{K_v}. It is necessary to set input "correction" in the function \code{\link[RKHSMetaMod]{calc_Kv}} equal to "TRUE".
}
\value{
List of \eqn{4} components: "mus", "qs", "mu", "res":
\item{mus}{Vector, values of the evaluated penalty parameters \eqn{\mu_g} in the RKHS group lasso algorithm until it achieves \eqn{\mu_{g}(q)}.}
\item{qs}{Vector, number of active groups associated with each value of \eqn{\mu_g} in mus.}
\item{mu}{Scalar, value of \eqn{\mu_{g}(q)}.}
\item{res}{An RKHS Group Lasso object:}
\item{intercept}{Scalar, estimated value of intercept.}
\item{teta}{Matrix with vMax rows and \eqn{n} columns. Each row of the matrix is the estimated vector \eqn{\theta_{v}} for \eqn{v=1,...,}vMax.}
\item{fit.v}{Matrix with \eqn{n} rows and vMax columns. Each row of the matrix is the estimated value of \eqn{f_{v}=K_{v}\theta_{v}}.}
\item{fitted}{Vector of size \eqn{n}, indicates the estimator of \eqn{m}.}
\item{Norm.H}{Vector of size vMax, estimated values of the penalty norm.}
\item{supp}{Vector of active groups.}
\item{Nsupp}{Vector of the names of the active groups.}
\item{SCR}{Scalar, equals to \eqn{\Vert Y-f_{0}-\sum_{v}K_{v}\theta_{v}\Vert ^{2}}.}
\item{crit}{Scalar, indicates the value of the penalized criteria.}
\item{MaxIter}{Integer, number of iterations until convergence is reached.}
\item{convergence}{TRUE or FALSE. Indicates whether the algorithm has converged or not.}
\item{RelDiffCrit}{Scalar, value of the first convergence criteria at the last iteration, \eqn{\frac{crit_{lastIter}-crit_{lastIter-1}}{crit_{lastIter-1}}}.}
\item{RelDiffPar}{Scalar, value of the second convergence criteria at the last iteration, \eqn{\Vert\frac{\theta_{lastIter}-\theta_{lastIter-1}}{\theta_{lastIter-1}}\Vert ^{2}}.}
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
\code{\link[RKHSMetaMod]{calc_Kv}}, \code{\link[RKHSMetaMod]{mu_max}}
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
result <- grplasso_q(Y,Kv,5,100 ,Num=10)
result$mu
result$res$Nsupp
}
\keyword{ Meta Model }
\keyword{RKHS Group Lasso}
\keyword{Penalty parameter }
\keyword{Active groups }
\keyword{ Optimization problem }