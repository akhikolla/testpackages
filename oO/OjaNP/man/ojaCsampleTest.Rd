\name{ojaCsampleTest}
\alias{ojaCsampleTest}
\alias{ojaCsampleTest.default}
\alias{ojaCsampleTest.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
C Sample Location Test Based on Oja Signs and Ranks
}
\description{
Function to test for equality of location in the C sample case using Oja signs and ranks.
}

\usage{ojaCsampleTest(X, \dots)

\method{ojaCsampleTest}{default}(X, Y, mu = NULL, scores = "sign", p = 1, 
               method = "approximation", n.simu = 1000, 
               center = "ojaMedian", na.action = na.fail, ...)

\method{ojaCsampleTest}{formula}(formula, scores="sign", p = 1, 
               method = "approximation", n.simu = 1000, 
               center = "ojaMedian", data, subset, na.action,...)
   }


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric data frame or matrix in the two sample case.
}
  \item{Y}{
a numeric data frame or matrix in the two sample case.
}
\item{formula}{a formula of the form \code{X ~ g} where \code{X}
    is a numeric matrix with at least two columns giving the data values and \code{g} a factor
    with at least two levels giving the corresponding groups.}
  \item{mu}{
a vector indicating the hypothesized value of the difference
in location. NULL represents no difference between the groups.
For more than two groups \code{mu} should be 0 or not be specified at all.
}
  \item{scores}{
options are \dQuote{rank} for the Oja rank test, \dQuote{sign} for the Oja sign test. The sign test is the default.
}
  \item{p}{
value of \dQuote{\code{p}} to be passed on to \code{\link{ojaSign}} or \code{\link{ojaRank}}. The default here is to use all hyperplanes since only then the tests are valid.
This can make the functions quite slow, however.
}
  \item{method}{
defines the method used for the computation of the p-value. The possibilities are
\dQuote{approximation} (default) or \dQuote{permutation}.
}
  \item{n.simu}{
if \dQuote{\code{method} = permutation} specifies this the number of replications used in the
                 permutation procedure.
}
  \item{center}{
value of \dQuote{\code{center}} to be passed on to \code{\link{ojaSign}}. Is used to center the data matrix. The default is the natural but computationally expansive
Oja median. For other options see the help for \code{\link{ojaSign}}.
}
  \item{data}{an optional data frame, list or environment  containing
          the variables in the model.  If not found in \dQuote{data}, the
          variables are taken from \dQuote{environment(formula)}.
  
  }
  \item{subset}{an optional vector specifying a subset of observations to be
          used for the testing.
  }
  \item{na.action}{
a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.
}
  \item{\dots}{
further arguments to be passed to or from methods. Mainly further arguments for \code{\link{ojaSign}} or \code{\link{ojaRank}}.
}
}
\details{
In the C-sample case of the Oja sign test the covariance matrix of the signs is divided by the sample size and not by sample size - 1.

For the sign test version always the Oja median should be used from a theoretical point of view to center the data and the median should be computed using
the exact algorithm. For further details about the the Oja median see \code{\link{ojaMedian}}.

Note that no theoretical results are available when \dQuote{\code{p}} is not set to 1.
}

\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the Q-statistic.}
  \item{parameter}{the degrees of freedom for the Q-statistic or the number of replications in the permutation procedure.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified hypothesized value of the difference in location.  (only in the two sample case)}
  \item{alternative}{a character string with the value 'two.sided'. (only in the two sample case)}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data.}
}
\references{
\cite{Hettmansperger, T. P. and Oja, H. (1994), Affine invariant multivariate multisample sign test, \emph{Journal of the Royal Statistical Society}, Series B, \bold{56}, 235--249.}

\cite{Hettmansperger, T. P., \enc{Möttönen}{Mottonen}, J. and Oja, H. (1999), Multivariate affine invariant rank tests for several samples, \emph{Statistica Sinica}, \bold{8}, 785--800.} 

\cite{Visuri, S., Ollila, E., Koivunen, V., \enc{Möttönen}{Mottonen}, J. and Oja, H. (2003), Affine equivariant multivariat rank methods, \emph{Journal of Statistical Planning and Inference}, \bold{114}, 161--185.} 

\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{ojaSign}}, \code{\link{ojaRank}}, \code{\link{oja1sampleTest}}
}
\examples{
data(biochem)
X <- subset(biochem, group=="Control", select=c("comp.1","comp.2"))
Y <- subset(biochem, group=="Treat", select=c("comp.1","comp.2"))
ojaCsampleTest(X,Y, alg="exact")
ojaCsampleTest(X,Y, method="p", alg="exact")
ojaCsampleTest(cbind(comp.1, comp.2) ~ group, score="r", data=biochem)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
