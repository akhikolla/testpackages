\name{mu_max}
\alias{mu_max}

\title{
Function to find the maximal value of the penalty parameter in the RKHS Group Lasso problem.
}
\description{
Calculates the value of the penalty parameter in the RKHS group lasso problem when the first penalized parameter group enters the model.
}
\usage{
mu_max(Y, matZ)
}
\arguments{
  \item{Y}{Vector of response observations of size \eqn{n}.}
  \item{matZ}{List of vMax components. Each component includes the eigenvalues and eigenvectors of the positive definite Gram matrices \eqn{K_v, v=1,...,}vMax. It should have the same format as the output "kv" of the function \code{\link[RKHSMetaMod]{calc_Kv}}.}
}
\details{
Details.
}
\value{
An object of type numeric is returned.
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
Kv <- calc_Kv(X, kernel, Dmax, TRUE,TRUE)
matZ <- Kv$kv
mumax <- mu_max(Y, matZ)
mumax
}
\keyword{RKHS Group Lasso }
\keyword{ Penalty parameter }
