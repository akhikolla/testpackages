\name{guo}
\alias{guo}
\title{
Mixed directional FDR controlled test decisions
}
\description{
Compute mixed directional FDR controlled test decisions for raw p-values
of multidimensional ordered tests
}
\usage{
guo(p,a)
}
\arguments{
  \item{p}{
    a numeric matrix, of raw p-values from multidimensional ordered
    tests. Each line corresponds to a group of ordered hypotheses
    e.g. t-tests between mean mesurements from consecutive levels of
    some ordered grouping variable. 
  }
  \item{a}{
    desired alpha level for the test
  }
}
\details{
guo implements the method proposed in Guo et al. (2010) that provides
test decisions controlling the mixed directional FDR for
multidimensional ordered tests. 
}
\value{
  A matrix with test decisions TRUE corresponds to significant rejection
  of the null and FALSE to acceptance of the null.
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

