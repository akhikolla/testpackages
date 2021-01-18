\name{pev_score}
\alias{pev_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PEV Score
}
\description{
PEV-score (Akdemir D. et al. (2015) <DOI:10.1186/s12711-015-0116-6>) is a criterion for finding a training set which derived from the covariance of the prediction of the test set.
}
\usage{
pev_score(X, X0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A genetic matrix of training set.
}
  \item{X0}{
A genetic matrix of test set.
}
}
%\details{}
\value{
A numeric score.
}
\references{
Akdemir D. et al., (2015), Optimization of genomic selection training populations with a genetic algoritthm. Genetics Selection Evolution 47:38.\\

Kennedy B.W., Trus D., (1993), Considerations on genetic connectedness between management units under an animal model. J Anim Sci. 1993 Sep;71(9):2341-52
}
\author{
Jen-Hsiang Ou <oumark.me@outlook.com>
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
\examples{
data("rice44kPCA")
pev_score(geno[1:50,],geno[51:100,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
