\name{r_score}
\alias{r_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
r Score
}
\description{
A criterion for finding training set which derived from Pearson's correlation between GEBVs (genomic estimated breeding value) and phenotype value of a test set.
}
\usage{
r_score(x, x0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A genetic matrix of training set.
}
  \item{x0}{
A genetic matrix of test set.
}
}
%\details{}
\value{
A numeric score.
}
\references{
Ou J.H., Liao C.T. (2018), Training set determination for genomic selection. National Taiwan University Master Thesis.
}
\author{
Jen-Hsiang Ou <oumark.me@outlook.com>
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
\examples{
data("rice44kPCA")
r_score(geno[1:50,],geno[51:100,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
