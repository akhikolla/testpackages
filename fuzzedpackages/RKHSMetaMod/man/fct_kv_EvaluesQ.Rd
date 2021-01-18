\name{calc_Kv}
\alias{calc_Kv}
\title{
Function to calculate the Gram matrices and their eigenvalues and eigenvectors for a chosen reproducing kernel.
}
\description{
Calculates the Gram matrices \eqn{K_v} for \eqn{v=1,\ldots,}vMax, and returns their associated eigenvalues and eigenvectors. The calculated Gram matrices may be not positive definite. The option "correction" of this function allows to replace the matrices \eqn{K_v} that are not positive definite by their "nearest positive definite" matrices.  
}
\usage{
calc_Kv(X, kernel, Dmax, correction, verbose, tol)
}
\arguments{
  \item{X}{Matrix of observations with \eqn{n} rows and \eqn{d} columns.}
  \item{kernel}{Character, the type of the reproducing kernel: matern \eqn{(}matern kernel\eqn{)}, brownian \eqn{(}brownian kernel\eqn{)}, gaussian \eqn{(}gaussian kernel\eqn{)}, linear \eqn{(}linear kernel\eqn{)}, quad \eqn{(}quadratic kernel\eqn{)}.}
  \item{Dmax}{Integer, between \eqn{1} and \eqn{d}, indicates the order of interactions considered in the meta model: Dmax\eqn{=1} is used to consider only the main effects, Dmax\eqn{=2} to include the main effects and the interactions of order \eqn{2,\ldots}.}
  \item{correction}{Logical, if TRUE, the program makes the correction to the matrices \eqn{K_v} that are not positive definite (see details). Set as TRUE by default.}
  \item{verbose}{Logical, if TRUE, the group \eqn{v} for which the correction is done is printed. Set as TRUE by default.}
  \item{tol}{Scalar, used if correction is TRUE. For each matrix \eqn{K_v} if \eqn{\lambda_{min} < \lambda_{max}\times}tol, then the correction to \eqn{K_v} is done (see details). Set as \eqn{1e^{-8}} by default.}
}
\details{
Let \eqn{\lambda_{v,i},i=1,...,n} be the eigenvalues associated with matrix \eqn{K_v}. Set \eqn{\lambda_{max}={max}_{i}\lambda_{v,i}} and \eqn{\lambda_{min}={min}_{i}\lambda_{v,i}}. The eigenvalues of \eqn{K_v} that is not positive definite are replaced by \eqn{\lambda_{v,i}+}epsilon, with espilon\eqn{=\lambda_{max}\times}tol. The value of tol depends on the type of the kernel and it is chosen small.
}
\value{
List of two components "names.Grp" and "kv":
\item{names.Grp}{Vector of size vMax, indicates the name of groups included in the meta model.}
\item{kv}{List of vMax components with the same names as the vector names.Grp. Each element of the list is a list of two components "Evalues" and "Q":}
\item{Evalues}{Vector of size \eqn{n}, eigenvalues of each Gram matrix \eqn{K_v}.}
\item{Q}{Matrix with \eqn{n} rows and \eqn{n} columns, eigenvectors of each Gram matrix \eqn{K_v}.}
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
Kv <- calc_Kv(X, kernel, Dmax)
names <- Kv$names.Grp
Eigen.val1 <- Kv$kv$v1.$Evalues
Eigen.vec1 <- Kv$kv$v1.$Q
}

\keyword{ Reproducing Kernel Hilbert Spaces }
\keyword{ Positive definite matrix }
