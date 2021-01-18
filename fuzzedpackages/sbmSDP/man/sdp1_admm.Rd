\name{sdp1_admm}
\alias{sdp1_admm}
\title{SDP-1 algorithm}
\description{Fits a balanced stochastic block model to an adjacency matrix using SDP-1. The function implements a first-order ADMM solver for SDP-1.}
\usage{sdp1_admm(As, K, opts)}
\arguments{
   \item{As}{ a binary adjacency matrix. }
   \item{K}{ number of communities (or blocks).}
   \item{opts}{ a list containing options. Pass the empty list, that is, "list()", to use the default values.  (See examples.)}
}
\value{
  A list containing the following items:
  \item{X}{ the estimated cluster matrix.}
  \item{delta}{ a vector of norm differences between consecutive cluster matrices at each step of the ADMM iteration.}
  \item{T_term}{ number of actual iterations performed.}
}
\references{
  On Semidefinite relaxations of the block model by A.A. Amini and E. Levina.
}
\author{
Arash A. Amini
}
\examples{
# Create a simple blkmodel with K=3 communities each of size m=20
blkmodel <- list(m=20, K=3, p=.9, q=.4)
blkmodel <- within(blkmodel, { 
                   n <- m*K
                   M <- kronecker(matrix(c(p,q,q,q,p,q,q,q,p),nrow=3),matrix(1,m,m))
                   As <- 1*(matrix(runif(n^2),nrow=n) < M)
                   })
# Call sdp1_admm with options:
#  rho  the ADMM parameter, 
#  T    maximum number of iteration
#  tol  tolerance for norm(X_{t+1} - X_t)
#  report_interval  how many iteration between reporting progress
sdp.fit <- with(blkmodel, 
                sdp1_admm(as.matrix(As), K, list(rho=.1, T=10000, tol=1e-5, report_interval=100)))

# plot the adjacency matrix and the estimated cluster matrix
par(mfrow=c(1,2))
image(blkmodel$As)
image(sdp.fit$X)
}