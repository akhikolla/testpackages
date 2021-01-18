\name{tPCAladle}
\alias{tPCAladle}

\title{
Ladle Estimate for tPCA
}
\description{
For r-dimensional tensors, the Ladle estimate for tPCA assumes that for a given mode \eqn{m}, the last \eqn{p_m - k_m} modewise eigenvalues are equal. Combining information from the eigenvalues and eigenvectors of the m-mode covariance matrix the ladle estimator yields estimates for \eqn{k_1,...,k_r}.}

\usage{
tPCAladle(x, n.boot = 200, ncomp = NULL)
}

\arguments{
  \item{x}{array of an order at least two with the last dimension corresponding to the sampling units.}
  \item{n.boot}{number of bootstrapping samples to be used.}
  \item{ncomp}{vector giving the number of components among which the ladle estimator is to be searched for each mode. The default follows the recommendation of Luo and Li 2016.}
}
\details{
The model here assumes that the eigenvalues of the m-mode covariance matrix are of the form \eqn{\lambda_{1,m} \geq ... \geq \lambda_{k_m,m} > \lambda_{k_m+1,m} =  ... = \lambda_{p_m,m}}{lambda_{1,m} >= ... >= lambda_{k_m,m} > lambda_{k_m+1,m} =  ... = \lambda_{p_m,m}}
and the goal is to estimate the value of \eqn{k_m} for all modes. The ladle estimate for this purpose combines the values of the 
scaled eigenvalues and the variation of the eigenvectors based on bootstrapping. The idea there is that for distinct eigenvales the variation of the eigenvectors
is small and for equal eigenvalues the corresponding eigenvectors have large variation.

This measure is then computed assuming \eqn{k_m}=0,..., \code{ncomp[m]} and the ladle estimate for \eqn{k_m} is the value where the measure takes its minimum. 

}


\value{
A list of class `tladle` containing:
  \item{U}{list containing the modewise rotation matrices.}
  \item{D}{list containing the modewise eigenvalues.}
  \item{S}{array of the same size as \code{x} containing the principal components.}
  \item{ResMode}{a list with the modewise results which are lists containing:
  \describe{
  \item{mode}{label for the mode.}
  \item{k}{the estimated value of k.}
  \item{fn}{vector giving the measures of variation of the eigenvectors using the bootstrapped eigenvectors for the different number of components.}
  \item{phin}{normalized eigenvalues.}
  \item{lambda}{the unnormalized eigenvalues used to compute phin.}
  \item{gn}{the main criterion for the ladle estimate - the sum of fn and phin. k is the value where gn takes its minimum.}
  \item{comp}{vector from 0 to the number of dimensions to be evaluated.}
  }}
  \item{xmu}{the data location}
  \item{data.name}{string with the name of the input data}
  \item{method}{string \code{tPCA}.}
  }


\references{
\cite{Koesner, C, Nordhausen, K. and Virta, J. (2019), Estimating the signal tensor dimension using tensorial PCA. Manuscript.}

\cite{Luo, W. and Li, B. (2016), Combining Eigenvalues and Variation of Eigenvectors for Order Determination, Biometrika, 103, 875--887. <doi:10.1093/biomet/asw051>}
}

\author{
Klaus Nordhausen
}



\seealso{
\code{\link{tPCA}}, \code{\link{ggtladleplot}}
}
\examples{
library(ICtest)
n <- 200
sig <- 0.6

Z <- rbind(sqrt(0.7)*rt(n,df=5)*sqrt(3/5),
           sqrt(0.3)*runif(n,-sqrt(3),sqrt(3)),
           sqrt(0.3)*(rchisq(n,df=3)-3)/sqrt(6),
           sqrt(0.9)*(rexp(n)-1),
           sqrt(0.1)*rlogis(n,0,sqrt(3)/pi),
           sqrt(0.5)*(rbeta(n,2,2)-0.5)*sqrt(20)
)

dim(Z) <- c(3, 2, n)

U1 <- rorth(12)[,1:3]
U2 <- rorth(8)[,1:2]
U <- list(U1=U1, U2=U2)
Y <- tensorTransform2(Z,U,1:2)
EPS <- array(rnorm(12*8*n, mean=0, sd=sig), dim=c(12,8,n))
X <- Y + EPS


TEST <- tPCAladle(X)
TEST
ggtladleplot(TEST)
}

\keyword{multivariate}
\keyword{array}