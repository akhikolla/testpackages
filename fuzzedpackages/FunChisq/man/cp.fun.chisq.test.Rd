\name{cp.fun.chisq.test}
\alias{cp.fun.chisq.test}
\title{
Comparative Chi-Squared Test for Model-Free Functional Heterogeneity
}
\description{
Comparative functional chi-squared tests on two or more contingency tables.
}

\usage{
cp.fun.chisq.test(
  x, method = c("fchisq", "nfchisq", "default", "normalized"),
  log.p = FALSE
)
}

\arguments{
  \item{x}{
  a list of at least two matrices representing contingency tables of the same dimensionality.
}

  \item{method}{
  a character string to specify the method to compute the functional chi-squared statistic and its p-value. The default is \code{"fchisq"} (equivalent to \code{"default"}). See Details.

    Note: \code{"default"} and \code{"normalized"} are deprecated.
}

  \item{log.p}{
  logical; if \code{TRUE}, the p-value is given as \code{log(p)}. Taking the log improves the accuracy when p-value is close to zero. The default is \code{FALSE}.
  }
}

\details{
The comparative functional chi-squared test determines whether the patterns underlying the contingency tables are heterogeneous in a functional way \insertCite{zhang2014nonparametric}{FunChisq}. Specifically, it evaluates whether the column variable is a changed function of the row variable across the contingency tables.

Two methods are provided to compute the functional chi-squared statistic and its p-value. When \code{method = "fchisq"} (or \code{"default"}), the p-value is computed using the chi-squared distribution; when \code{method =} \code{"nfchisq"} (or \code{"normalized"}) a normalized statistic is obtained by shifting and scaling the original statistic and a p-value is computed using the standard normal distribution \insertCite{Box2005}{FunChisq}
 (Box et al., 2005). The normalized test is more conservative on the degrees of freedom.
}

\value{
A list with class "\code{htest}" containing the following components:

\item{statistic}{functional heterogeneity statistic if \code{method = "fchisq"} (equivalent to \code{"default"}), or normalized statistic if \code{method = "nfchisq"} (equivalent to \code{"normalized"}).}
\item{parameter}{degrees of freedom.}
\item{p.value}{p-value of the comparative functional chi-squared test. By default, it is computed by the chi-squared distribution. If \code{method = "normalized"}, it is the p-value of the normalized statistic computed by the standard normal distribution.}
}

\references{
\insertAllCited{}
}

\author{
Yang Zhang and Joe Song
}

\seealso{
For comparative chi-squared test that does not consider functional dependencies, \code{\link[DiffXTables]{cp.chisq.test}}.
}

\examples{
\donttest{
x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
y <- t(x)
z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
data <- list(x,y,z)
cp.fun.chisq.test(data)
cp.fun.chisq.test(data, method="nfchisq")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{htest}
\keyword{nonparametric}
