\name{e2test}
\alias{e2test}
\title{
Perform permutation Barlow's test on a dataset
}
\description{
This function computes Barlow's test for each line of a given
matrix. The global null distribution is computed using permutation. FWE
control is provided by the maxT procedure.
}
\usage{
e2test(data,g,B,rep=rep(1,length(g)))
}
\arguments{
  \item{data}{
    a numeric matrix, for the lines of which we want to test the null of no
    trend across groups
  }
  \item{g}{
    integer vector with group labels specified as ordered integers
    ranging from 1 to n where n is the number of ordered categories
  }
  \item{B}{
    number of permutations or a permutation matrix
  }
  \item{rep}{
    integer vector with group labels ranging from 1 to m where m is the
    number of independent samples (e.g. individuals). The ordering is
    not important.
  }
}
\details{
e2test takes a matrix for each line of which the null hypotheses of no
trend across ordered groups is tested using a permutation test based on
Barlow's E2 statistic. By permuting only samples within independent
entities of the experimental design, random factors such as technical
replicates from the same sample material can be incorporated. Multiple
testing control is provided by the maxT procedure. 
}
\value{
  \item{adj}{
    vector with FWE adjusted p.values for each null hypothesis}
  \item{raw}{
    vector with raw p.values}
  \item{dir}{
    vector with the direction of each decision
    }
}
\references{
  Klinglmueller, F., Tuechler, T., Posch, M. (2010)
  "Cross Platform Comparison Of Microarray Data Using Order
    Restricted Inference" Under Review
  
  Barlow, R. E., Bartholomew, D. J., Bremner, J. M., and Brunk, H. D. (1972)
  "Statistical inference under order restrictions"; Wiley, London.

  Robertson, T., Wright,F. T. and Dykstra, R. L. (1988)
  "Order Restricted Statistical Inference"; Wiley, New York.
}
\author{
  Florian Klinglmueller <float_at_lefant.net>
  Part of this function is C code that has been ported by Korbinian Strimmer from
  R code originally written by Kaspar Rufibach. 
}
\examples{
data <- matrix(rnorm(7200),nc=72)
groups <- rep(1:4,each=18)
ind <- rep(rep(1:3,each=6),4)
out <- e2test(data,groups,B=1000,rep=ind)
sum(out$adj<.05)
data2 <- data+matrix(rep(groups,nrow(data)),nr=nrow(data),byrow=TRUE)
out2 <- e2test(data2,groups,B=1000,rep=ind)
sum(out2$adj<.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }

