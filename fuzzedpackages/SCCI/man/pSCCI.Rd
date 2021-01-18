\name{pSCCI}
\alias{pSCCI}
\docType{package}
\title{
Stochastic Complexity-based Conditional Independence Criterium (p-value)
}
\description{
This is an adapted version of \eqn{SCCI} for which the output can be interpreted as a p-value. For this, we adapted \eqn{SCCI} such that if \eqn{SCCI = 0} (\eqn{X} is independent of \eqn{Y} given \eqn{Z}) it gives a p-value greater than \eqn{0.01} and for \eqn{SCCI > 0} (\eqn{X} is not independent of \eqn{Y} given \eqn{Z}) gives a p-value smaller or equal to \eqn{0.01}. Note that we just transformed the output of \eqn{SCCI} and do not obtain a real p-value. In essence, we define the artificial p-value as follows. Let \emph{v} the output of \eqn{SCCI} divided by the number of samples \emph{n}. \eqn{p = 2^{-(6.643855 - v)}}, which is equal to \eqn{0.01000001} if \eqn{v = 0}. Further, \eqn{p \le 0.01} for \eqn{SCCI \ge 0.000001}. We restrict the p-values to be between \eqn{0} and \eqn{1}.

Unlike \eqn{SCCI}, \eqn{pSCCI} is currently only instantiated with fNML.

\eqn{pSCCI} can be used directly in the PC algorithm developed by Spirtes et al. (2000), which was implemented in the 'pcalg' R-package by Kalisch et al. (2012), as shown in the example.
}
\arguments{
  \item{x}{Position of x variabe (integer).}
  \item{y}{Position of y variabe (integer).} 
  \item{S}{Vector of the position of zero or more conditioning variables (integer).}
  \item{suffStat}{This format was adapted such that it can be used in the PC algorithm and other algorithms from the 'pcalg' package. \eqn{SCCI} only need the filed "dm" that contains the data matrix.}
}
\usage{
pSCCI(x, y, S, suffStat)	
}
\examples{
set.seed(1)
x = round((runif(1000, min=0, max=5)))
y = round((runif(1000, min=0, max=5)))
Z = data.frame(round((runif(1000, min=0, max=5))), round((runif(1000, min=0, max=5))))
## create data matrix
data_matrix = as.matrix(data.frame(x,y,S1=Z[,1], S2=Z[,2]))
suffStat = list(dm=data_matrix)
pSCCI(x=1,y=2,S=c(3,4),suffStat=suffStat)	## 0.01000001

### Using SCI within the PC algorithm
if(require(pcalg)){
  ## Load data
  data(gmD)
  V <- colnames(gmD$x)
  ## define sufficient statistics
  suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)
  ## estimate CPDAG
  pc.D <- pc(suffStat,
            ## independence test: SCCI using fNML
            indepTest = pSCCI, alpha = 0.01, labels = V, verbose = TRUE)
}
if (require(pcalg) & require(Rgraphviz)) {
  ## show estimated CPDAG
  par(mfrow = c(1,2))
  plot(pc.D, main = "Estimated CPDAG")
  plot(gmD$g, main = "True DAG")
}
}
\references{
  Markus Kalisch, Martin Mächler, Diego Colombo, Marloes H. Maathuis, Peter Bühlmann; Causal inference using graphical models with the R package pcalg, Journal of Statistical Software, 2012

  Alexander Marx and Jilles Vreeken; Testing Conditional Independence on Discrete Data using Stochastic Complexity, Proceedings of the 22nd International Conference on Artificial Intelligence and Statistics (AISTATS), PMLR, 2019
  
  Peter Spirtes, Clark N. Glymour, Richard Scheines, David Heckerman, Christopher Meek, Gregory Cooper and Thomas Richardson; Causation, Prediction, and Search, MIT press, 2000
}