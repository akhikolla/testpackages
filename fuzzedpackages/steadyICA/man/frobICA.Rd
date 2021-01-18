\name{frobICA}
\alias{frobICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{match mixing matrices or ICs and calculate their Frobenius distance}
\description{
The ICA model is only identifiable up to signed permutations of the ICs. This function provides a similarity measure between two mixing matrices for the model X = S M + E, where X is n x p, S is n x d, and M is d x p. The input is either two mixing matrices M1 and M2 or two matrices of independent components S1 and S2. For M1 and M2, frobICA() finds the signed row permutation of M2 that minimizes the Frobenius norm between M1 and M2 using the Hungarian method. For S1 and S2, frobICA() finds the signed column permutation of S2 that minimizes the Frobenius norm between S1 and S2. This function allows the mixing matrices (or independent components) to have differing numbers of rows (respectively, columns) such that the similarity measure is defined by the matching rows (resp., columns), and the non-matching rows (resp., columns) are discarded.}
\usage{
frobICA(M1 = NULL, M2 = NULL, S1 = NULL, S2 = NULL, standardize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M1}{A d x p mixing matrix}
  \item{M2}{A d x q mixing matrix}
  \item{S1}{An n x d matrix of independent components}
  \item{S2}{An n x q matrix of independent components}
  \item{standardize}{Logical. See Note.}
}

\details{
frobICA(M1,M2) = 0 if there exists a signed permutation of the rows of M2 such that M1 = P\%*\%M2, where P is a d x q signed permutation matrix, i.e., composed of 0, 1, and -1, with d <= q; the function also allows d > q, in which case frobICA(M1,M2) = 0 if there exists a P such that P\%*\% M1 = M2. Unlike other ICA performance measures, this function can accomodate non-square mixing matrices.
}
\value{returns the Frobenius norm divided by p*min(d,q) (or n*min(d,q)) of the matched mixing matrices (resp., matched independent components).}
\note{
 If standardize=TRUE, then scales the rows of M1 and M2 to have unit norm or the columns of S1 and S2 to have zero mean and sample variance equal to one. The user can supply either M1 and M2 or S1 and S2 but not both.}
\references{
Kuhn, H. The Hungarian Method for the assignment problem Naval Research Logistics Quarterly, 1955, 2, 83 - 97

Risk, B.B., D.S. Matteson, D. Ruppert, A. Eloyan, B.S. Caffo. In review, 2013. Evaluating ICA methods with an application to resting state fMRI. 
}
\author{
Benjamin Risk
}

\seealso{
\code{\link[JADE:MD]{JADE::MD}}
\code{\link[clue:solve_LSAP]{clue::solve_LSAP}}
\code{\link{matchICA}}
}

\examples{
mat1 <- matrix(rnorm(4*6),nrow=4)
perm <- matrix(c(-1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1),4,4)
mat2 <- perm\%*\%mat1
sqrt(sum((mat1-mat2)^2))
frobICA(M1=mat1,M2=mat2)

#Another example showing invariance to permutations:
covMat <- t(mat1)\%*\%mat1
mvsample <- matrix(rnorm(400),100,4)\%*\%mat1
frobICA(M1=cov(mvsample),M2=covMat)
frobICA(M1=cov(mvsample),M2=covMat[sample(1:6),])

#Example using independent components:
nObs=300
simS<-cbind(rgamma(nObs, shape = 1, scale = 2),
            rgamma(nObs, shape = 3, scale = 2),
            rgamma(nObs, shape = 3, scale = 2),
            rgamma(nObs, shape = 9, scale = 0.5))
              
#not necessary in this example, but this should be done when used with ICA:            
simS <- apply(simS,2,scale) 
frobICA(S1=simS,S2=simS\%*\%perm)
\dontrun{
#returns an error if S1 and S2 are not explicitly defined:
frobICA(simS,simS\%*\%perm)            
}            

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ICA}
\keyword{distance}% __ONLY ONE__ keyword per line
