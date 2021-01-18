\name{pttest}
\alias{pttest}
\title{
Perform permutation t-tests between consecutive levels of ordered groups
}
\description{
This function computes permutation t-tests for differences in each line
of a matrix between consecutive levels of an ordered grouping variable.
}
\usage{
pttest(data,g,B,rep=rep(1,length(g)))
}
\arguments{
  \item{data}{
    a numeric matrix, for the lines of which we want to test the null of no
    trend across groups
  }
  \item{g}{
    integer vector with group labels for the ordered categories coded as
    integers between 1 and n, where n is the number of ordered categories
  }
  \item{B}{
    number of permutations or a list of permutation matrices
  }
  \item{rep}{
    integer vector with group labels ranging from 1 to m where m is the
    number of independent samples (e.g. individuals). The ordering is
    not important.
  }
}
\details{
pttest takes a matrix and for each line computes a permutation t-test
between consecutive levels of an ordered grouping variable. Random
factors e.g. technical replicates can be specified and will be accounted
for by permuting samples only within independent units of this
factor. A matrix with raw permutation p-values as well as a matrix with
the directional decisions are returned. See \code{\link{guo}} and
\code{\link{sigdirPttest}}for a way to get 
test decisions with control over the mixed directional false discovery
rate. 
}
\value{
  \item{raw}{
    matrix with raw p.values}
  \item{dir}{
    matrix with the direction of each decision
    }
}
\references{
  Klinglmueller, F., Tuechler, T., Posch, M. (2010)
  "Cross Platform Comparison Of Microarray Data Using Order
    Restricted Inference" Under Review

  Guo W., Sarkar SK., Peddada SD. (2010)
  "Controlling False Discoveries in Multidimensional Directional
  Decisions, with Applications to Gene Expression Data on Ordered
  Categories" Biometrics
}
\author{
  Florian Klinglmueller <float_at_lefant.net>
}
\examples{
data <- matrix(rnorm(7200),nc=72)
groups <- rep(1:4,each=18)
ind <- rep(rep(1:3,each=6),4)
out <- pttest(data,groups,B=1000,rep=ind)
guo(out$unadj,.05)
data2 <- data+matrix(rep(groups,nrow(data)),nr=nrow(data),byrow=TRUE)
out2 <- pttest(data2,groups,B=1000,rep=ind)
guo(out2$unadj,.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

