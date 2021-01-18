\name{RKHSMetMod_qmax}
\alias{RKHSMetMod_qmax}

\title{
Function to produce a sequence of meta models, with at most qmax active groups in each meta model. The meta models are the solutions of the RKHS Ridge Group Sparse or RKHS Group Lasso optimization problems.
}
\description{
Calculates the Gram matrices \eqn{K_v} for a chosen kernel, determines \eqn{\mu}, note \eqn{\mu (qmax)}, for which the number of active groups in the RKHS group lasso solution is equal to qmax, and fits a solution of an RKHS ridge group sparse or an RKHS group lasso problem for each pair of penalty parameters \eqn{(\mu (qmax),\gamma)}, in the Gaussian regression model.
}
\usage{
RKHSMetMod_qmax(Y, X, kernel, Dmax, gamma, qmax, rat, Num, verbose)
}
\arguments{
  \item{Y}{Vector of response observations of size \eqn{n}.}
  \item{X}{Matrix of observations with \eqn{n} rows and \eqn{d} columns.}
  \item{kernel}{Character, indicates the type of the reproducing kernel: matern \eqn{(}matern kernel\eqn{)}, brownian \eqn{(}brownian kernel\eqn{)}, gaussian \eqn{(}gaussian kernel\eqn{)}, linear \eqn{(}linear kernel\eqn{)}, quad \eqn{(}quadratic kernel\eqn{)}. See the function \code{\link[RKHSMetaMod]{calc_Kv}}.}
  \item{Dmax}{Integer, between \eqn{1} and \eqn{d}, indicates the order of interactions considered in the meta model: Dmax\eqn{=1} is used to consider only the main effects, Dmax\eqn{=2} to include the main effects and the interactions of order \eqn{2,\ldots}. See the function \code{\link[RKHSMetaMod]{calc_Kv}}.}
  \item{gamma}{Vector of non negative scalars, values of the penalty parameter \eqn{\gamma} in decreasing order. If \eqn{\gamma=0} the function solves an RKHS Group Lasso problem and for \eqn{\gamma>0} it solves an RKHS Ridge Group Sparse problem.}
  \item{qmax}{Integer, shows the maximum number of active groups in the obtained solution.}
  \item{rat}{Positive scalar, to restrict the minimum value of \eqn{\mu} considered in the algorithm, \eqn{\mu_{min}=\mu_{max}/(\sqrt{n}\times rat)}. The value \eqn{\mu_{max}} is calculated inside the program, see function \code{\link[RKHSMetaMod]{mu_max}}.}
  \item{Num}{Integer, it is used to restrict the number of different values of the penalty parameter \eqn{\mu} to be evaluated in the RKHS Group Lasso algorithm until it achieves \eqn{\mu (qmax)}: for Num \eqn{= 1} the program is done for \eqn{3} different values of \eqn{\mu}, \eqn{\mu_{1}=(\mu_{min}+\mu_{max})/2}, \eqn{\mu_{2}=(\mu_{min}+\mu_{1})/2} or \eqn{\mu_{2}=(\mu_{1}+\mu_{max})/2} depending on the number of active groups in the meta model associated with \eqn{\mu_{1}}, \eqn{\mu_{3}=\mu_{min}}.}
  \item{verbose}{Logical, if TRUE, prints: the group \eqn{v} for which the correction of Gram matrix \eqn{K_v} is done, and for each pair of \eqn{(\mu,\gamma)}: the number of current iteration, active groups and convergence criterias. Set as FALSE by default.}
}
\details{
Details.
}
\value{
List of three components "mus", "qs", and "MetaModel":
\item{mus}{Vector, values of the evaluated penalty parameters \eqn{\mu} in the RKHS Group Lasso algorithm until it achieves \eqn{\mu (qmax)}.}
\item{qs}{Vector, number of active groups associated with each element in mus.}
\item{MetaModel}{List with the same length as the vector gamma. Each component of the list is a list of \eqn{3} components "mu", "gamma" and "Meta-Model": }
\item{mu}{Scalar, the value \eqn{\mu (qmax)}.}
\item{gamma}{Positive scalar, element of the input vector gamma associated with the estimated Meta-Model.}
\item{Meta-Model}{An RKHS Ridge Group Sparse or RKHS Group Lasso object associated with the penalty parameters mu and gamma:}
\item{intercept}{Scalar, estimated value of intercept.}
\item{teta}{Matrix with vMax rows and \eqn{n} columns. Each row of the matrix is the estimated vector \eqn{\theta_{v}} for \eqn{v=1,...,}vMax.}
\item{fit.v}{Matrix with \eqn{n} rows and vMax columns. Each row of the matrix is the estimated value of \eqn{f_{v}=K_{v}\theta_{v}}.}
\item{fitted}{Vector of size \eqn{n}, indicates the estimator of \eqn{m}.}
\item{Norm.n}{Vector of size vMax, estimated values for the Ridge penalty norm.}
\item{Norm.H}{Vector of size vMax, estimated values of the Sparse Group penalty norm.}
\item{supp}{Vector of active groups.}
\item{Nsupp}{Vector of the names of the active groups.}
\item{SCR}{Scalar, equals to \eqn{\Vert Y-f_{0}-\sum_{v}K_{v}\theta_{v}\Vert ^{2}}.}
\item{crit}{Scalar, indicates the value of penalized criteria.}
\item{gamma.v}{Vector, coefficients of the Ridge penalty norm, \eqn{\sqrt{n}\gamma\times}gama_v.}
\item{mu.v}{Vector, coefficients of the Group Sparse penalty norm, \eqn{n\mu\times}mu_v.}
\item{iter}{List of two components: maxIter, and the number of iterations until the convergence is achieved.}
\item{convergence}{TRUE or FALSE. Indicates whether the algorithm has converged or not.}
\item{RelDiffCrit}{Scalar, value of the first convergence criteria at the last iteration, \eqn{\Vert\frac{\theta_{lastIter}-\theta_{lastIter-1}}{\theta_{lastIter-1}}\Vert ^{2}}.}
\item{RelDiffPar}{Scalar, value of the second convergence criteria at the last iteration, \eqn{\frac{crit_{lastIter}-crit_{lastIter-1}}{crit_{lastIter-1}}}.}
}
\references{
Kamari, H., Huet, S. and Taupin, M.-L. (2019) RKHSMetaMod : An R package to estimate the Hoeffding decomposition of an unknown function by solving RKHS Ridge Group Sparse optimization problem. <arXiv:1905.13695>
}
\author{
Halaleh Kamari
}
\note{
For the case \eqn{\gamma=0} the outputs "mu"\eqn{=\mu_{g}} and "Meta-Model" is the same as the one returned by the function \code{\link[RKHSMetaMod]{RKHSgrplasso}}.
}
\seealso{
\code{\link[RKHSMetaMod]{calc_Kv}}, \code{\link[RKHSMetaMod]{mu_max}}, \code{\link[RKHSMetaMod]{RKHSgrplasso}}, \code{\link[RKHSMetaMod]{pen_MetMod}}, \code{\link[RKHSMetaMod]{grplasso_q}}
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
gamma <- c(.5,.01,.001,0)
Num <- 10
rat <- 100
qmax <- 4
result <- RKHSMetMod_qmax(Y, X, kernel, Dmax, gamma, qmax, rat, Num,FALSE)
names(result)
result$mus
result$qs
l <- length(gamma)
for(i in 1:l){print(result$MetaModel[[i]]$mu)}
for(i in 1:l){print(result$MetaModel[[i]]$gamma)}
for(i in 1:l){print(result$MetaModel[[i]]$`Meta-Model`$Nsupp)}
}
\keyword{ Meta Model }
\keyword{ RKHS Ridge Group Sparse}
\keyword{ RKHS Group Lasso}
\keyword{ Optimization problem }
\keyword{ Reproducing Kernel Hilbert Spaces }
\keyword{Active groups }
\keyword{Penalty parameter }