\name{SCCI}
\alias{SCCI}
\docType{package}
\title{
Stochastic Complexity-based Conditional Independence Criterium
}
\description{
Calculates whether two random variables \eqn{X} and \eqn{Y} are independent given a set of variables \eqn{Z} using \eqn{SCCI}. A score of \eqn{0} denotes that independence holds and values greater than \eqn{0} mean that \eqn{X} is not independent of \eqn{Y} given \eqn{Z}. For details on \eqn{SCCI}, we refer to Marx and Vreeken (AISTATS, 2019). If you use \eqn{SCCI} in your work, please cite Marx and Vreeken (AISTATS, 19).

The output of \eqn{SCCI(..)} is the difference in \emph{number of bits} between condtioning \eqn{X} only on \eqn{Z} and conditioning on \eqn{Z} and \eqn{Y}. For the variant of \eqn{SCCI} that gives outputs that can be intpreted as p-values, please refer to pSCCI. 
}
\arguments{
  \item{x}{A discrete vector.}
  \item{y}{A discrete vector.} 
  \item{Z}{A data frame consisting of zero or more columns of discrete vectors.}
  \item{score}{Default: fNML, optionally qNML can be passed.}
  \item{sym}{sym can be true or false}
}
\usage{
SCCI(x, y, Z, score="fNML", sym=FALSE)	
}
\examples{
set.seed(1)
x = round((runif(1000, min=0, max=5)))
y = round((runif(1000, min=0, max=5)))
Z = data.frame(round((runif(1000, min=0, max=5))), round((runif(1000, min=0, max=5))))
SCCI(x=x,y=y,Z=Z,score="fNML",sym=FALSE)	## 0
}
\references{
  Alexander Marx and Jilles Vreeken; Testing Conditional Independence on Discrete Data using Stochastic Complexity, Proceedings of the 22nd International Conference on Artificial Intelligence and Statistics (AISTATS), PMLR, 2019
}