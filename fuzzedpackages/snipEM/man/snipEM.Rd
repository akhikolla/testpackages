\name{snipEM}
\alias{snipEM}
\alias{snipEM.initialV}
\title{Snipping for location and scatter estimation with cellwise outliers}
\description{
  Computes an estimator optimizing the Gaussian likelihood over a
  snipping set. The function \code{snipEM.initialV} can be used to
  perform some iterations to initialize \code{V}. 
}
\usage{
snipEM(X, V, tol = 1e-04, maxiters = 500, maxiters.S = 1000, print.it = FALSE) 

snipEM.initialV(X, V, mu0, S0, maxiters.S = 100, greedy = TRUE)
}
\arguments{
  \item{X}{Data.}
  \item{V}{Binary matrix of the same size as \code{X}. Zeros correspond to initial snipped entries.}
  \item{tol}{Tolerance for convergence. Default is \code{1e-4}.}
  \item{maxiters}{Maximum number of iterations for the SM algorithm. Default is \code{500}.}
  \item{maxiters.S}{Maximum number of iterations of the inner greedy snipping algorithm. Default is \code{1000}.}
  \item{print.it}{Logical; if \code{TRUE}, partial results are print. Default is \code{FALSE}.}
  \item{mu0}{Initial estimate for the mean vector that is used in the initialization stage.}
  \item{S0}{Initial estimate for the covariance matrix that is used in the initialization stage.}
  \item{greedy}{Logical; if \code{TRUE}, perform the greedy snipping algorithm in search for the binary 
  matrix that gives the largest likelihood value throughout \code{maxiters.S} iterations. 
  If \code{FALSE}, stop right after the snipping algorithm finds a binary matrix that gives a larger 
  likelihood value than the initial one. Default is \code{TRUE}.}

}
\details{
   This function computes the \code{sclust} estimator of Farcomeni
   (2014) with \eqn{k=1}. It therefore provides a robust estimate of
   location and scatter in presence of entry-wise outliers. It is
   based on a snip-maximize (SM) algorithm. At the S step, the
   likelihood is optimized over the set of snipped entries, at the M
   step the location and scatter estimates are updated. The S step is
   based on a greedy algorithm, unlike the one proposed in Farcomeni
   (2014,2014a). The number of snipped entries \code{sum(1-V)} is kept
   fixed throughout. 

   Results depend on good initialization of the \code{V} matrix. A
   boxplot rule (see examples) usually works well. The function
   \code{snipEM.initialV} can be used to improve the initial choice
   through some iterations updating only \code{V} from initial
   (robust) estimates \code{mu0} and \code{S0}. In the example, the
   EMVE is used to obtain \code{mu0} and \code{S0}.
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{mu} \tab Estimated location. \cr
	\code{S} \tab Estimated scatter matrix. \cr 
	\code{V} \tab Final (optimal) V matrix. \cr 
	\code{lik} \tab Gaussian log-likelihood at convergence. \cr 
	\code{iter} \tab Number of outer iterations before convergence. \cr 
	}
}
\references{
Farcomeni, A. (2014) Snipping for robust k-means clustering under component-wise contamination, \emph{Statistics and Computing}, \bold{24}, 909-917

Farcomeni, A. (2014) Robust constrained clustering in presence of entry-wise outliers, \emph{Technometrics}, \bold{56}, 102-111   
}
\author{Alessio Farcomeni \email{alessio.farcomeni@uniroma1.it}, Andy
Leung \email{andy.leung@stat.ubc.ca}}
\seealso{
  \code{\link{sclust}}, \code{\link{stEM}},
    \code{\link{sumlog}},
  \code{\link{ldmvnorm}} 
}
\examples{
n=100
p=5
Xc <- matrix(rnorm(100*10),100,5)

# initial V 
V <- matrix(1,n,p)
V[!is.na(match(as.vector(Xc),boxplot(as.vector(Xc),plot=FALSE)$out))] <- 0
Xna <- Xc 
Xna[ which( V == 0) ] <- NA

resSEM <- snipEM(Xc, V)

}


