\name{rightskew}
\alias{rightskew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
force ICs to have positive skewness and order by skewness
}
\description{
The ICA model is only identifiable up to signed permutations. This function provides a canonical ordering for ICA that is useful for fMRI or studies where signals are skewed. Multiplies columns of S that are left-skewed by -1 to force right skewness. Optionally orders the columns by descending skewness.
}
\usage{
rightskew(S, M = NULL, order.skew = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
n x d matrix 
}
  \item{M}{
d x p mixing matrix
}
  \item{order.skew}{
Option to return the permutation of columns of S from largest to smallest skewness. Also returns a permuted version of M that corresponds with the permuted S.
}
}
\value{
Returns the matrix S such that all columns have positive skewness. If optional argument M is supplied, returns a list with the new S and corresponding M.
}
\references{
Eloyan, A. & Ghosh, S. A Semiparametric Approach to Source Separation using Independent Component Analysis Computational Statistics and Data Analysis, 2013, 58, 383 - 396.
}
\author{
Benjamin Risk
}

\examples{
nObs = 1024
simS<-cbind(rgamma(nObs, shape = 1, scale = 2),
            rgamma(nObs, shape = 9, scale = 0.5),
            -1*rgamma(nObs, shape = 3, scale = 2))

apply(simS,2,function(x){ 
  (sum((x - mean(x))^3)/length(x))/(sum((x - mean(x))^2)/length(x))^(3/2)})

canonicalS <- rightskew(simS)

apply(canonicalS,2,function(x){
 (sum((x - mean(x))^3)/length(x))/(sum((x - mean(x))^2)/length(x))^(3/2)})

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
