\name{oja1sampleTest}
\Rdversion{1.1}
\alias{oja1sampleTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  One Sample Location Test Based on Oja Signs and Ranks
}
\description{
 Function to test for location in the one sample case using Oja signs and ranks. 
}
\usage{
oja1sampleTest(X, mu = NULL, scores = "sign", p = 1,
               method = "approximation", n.simu = 1000, 
               na.action = na.fail, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a numeric data frame or data matrix.}
  \item{mu}{a vector indicating the hypothesized value of the location. NULL represents the origin.}
  \item{scores}{options are \dQuote{rank} for the signed rank test, \dQuote{sign} for the sign test. The sign test is the default.}
  \item{p}{value of \dQuote{p} to be passed on to \code{\link{ojaSign}} or \code{\link{ojaSignedRank}}. The default here is to use all hyperplanes since only then the tests are valid.
           This can make the functions quite slow, however.} 
  \item{method}{defines the method used for the computation of the p-value. The possibilities are
        \dQuote{approximation} (default) or \dQuote{permutation} which permutes the signs of the Oja signs or Oja signed ranks.}
  \item{n.simu}{if \dQuote{\code{method} = permutation} specifies this the number of replications used in the
                 permutation procedure.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
  \item{\dots}{further arguments to be passed to or from methods.}
}

\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the Q-statistic.}
  \item{parameter}{the degrees of freedom for the Q-statistic or the number of replications in the permutation procedure.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified hypothesized value of the location.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data.}
}
\references{
\cite{Hettmansperger, T. P., Nyblom, J. and Oja, H. (1994), Affine invariant multivariate one-sample sign test, \emph{Journal of the Royal Statistical Society}, Series B, \bold{56}, 221--234.}

\cite{Hettmansperger, T. P., \enc{Möttönen}{Mottonen}, J. and Oja, H. (1997), Multivariate affine invariant one-sample signed-rank tests, \emph{Journal of the American Statistical Society}, \bold{92}, 1591--1600.} 

\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{ojaMedian}}, \code{\link{ojaSign}}, \code{\link{ojaSignedRank}}
}
\examples{
data(biochem)
oja1sampleTest(biochem[,1:2], mu = c(1.1, 0.4))
oja1sampleTest(biochem[,1:2], mu = c(1.1, 0.4), method = "p")

oja1sampleTest(biochem[,1:2], mu = c(1.1, 0.4), scores = "rank")
oja1sampleTest(biochem[,1:2], mu = c(1.1, 0.4), scores = "rank", method = "p")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
