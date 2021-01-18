\name{cd_score}
\alias{cd_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized Coefficient of Determination
}
\description{
A criterion for finding optimal training set using generalized coefficient of determination (Laloe D. (1993) <DOI:10.1186/1297-9686-25-6-557>).
}
\usage{
cd_score(x, x0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Genomic metrix of training set.
}
  \item{x0}{
Genomic metric of testing set.
}
}
%\details{}
\value{
A numeric score.
}
\references{
Laloe D. (1993), Precision and information in linear models of genetic evaluation. Genetics Selection Evolution 25:557.
}
\author{
Jen-Hsiang Ou <oumark.me@outlook.com>
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}
\examples{
data("rice44kPCA")
cd_score(geno[1:100,], geno[101:200,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
