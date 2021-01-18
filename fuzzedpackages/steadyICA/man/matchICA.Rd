\name{matchICA}
\alias{matchICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
match independent components using the Hungarian method
}
\description{
The ICA model is only identifiable up to signed permutations of the ICs. This function finds the signed permutation of a matrix S such that ||S\%*\%P - template|| is minimized. Optionally also matches the mixing matrix M.
}
\usage{
matchICA(S, template, M = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
the n x d matrix of ICs to be matched
}
  \item{template}{
the n x d matrix that S is matched to.
}
  \item{M}{
an optional d x p mixing matrix corresponding to S that will also be matched to the template
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns the signed permutation of S that is matched to the template. If the optional argument M is provided, returns a list with the permuted S and M matrices. 
}
\references{
Kuhn, H. The Hungarian Method for the assignment problem Naval Research Logistics Quarterly, 1955, 2, 83 - 97

Risk, B.B., D.S. Matteson, D. Ruppert, A. Eloyan, B.S. Caffo. In review, 2013. Evaluating ICA methods with an application to resting state fMRI. 
}
\author{
Benjamin Risk
}

\seealso{
\code{\link{frobICA}}
\code{\link[clue:solve_LSAP]{clue::solve_LSAP}}
}
\examples{
set.seed(999)
nObs <- 1024
nComp <- 3
# simulate from some gamma distributions:
simS<-cbind(rgamma(nObs, shape = 1, scale = 2),
            rgamma(nObs, shape = 3, scale = 2),
            rgamma(nObs, shape = 9, scale = 0.5))

simM <- matrix(rnorm(9),3)
pMat <- matrix(c(0,-1,0,1,0,0,0,0,-1),3)
permS <- simS\%*\%pMat
permM <- t(pMat)\%*\%simM

matchedS <- matchICA(S = permS, template = simS, M = permM)
sum(abs(matchedS$S - simS))
sum(abs(simM - matchedS$M))
}
\keyword{ hungarian }
\keyword{ ICA }% __ONLY ONE__ keyword per line
