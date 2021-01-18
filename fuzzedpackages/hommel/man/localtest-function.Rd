\name{Local test}

\alias{localtest}

\title{P-values and adjusted p-values of local tests}

\description{An intersection hypothesis is defined for any collection of null hypotheses as the hypothesis that is true if and only if the proportion of null hypotheses are less than the selected threshold. A test for an intersection hypothesis is called a local test. This function calculates familywise error adjusted p-values for local tests. Familywise error control is over all possible local tests simultaneously.}

\usage{
  localtest (hommel, ix, tdp)
}

\arguments{
  \item{hommel}{The hommel object for the complete multiple testing problem.}
  \item{ix}{The selection of hypotheses defining the intersection hypothesis. Any selector that works on the original vector of p-values (index, negative index, logical, names) is allowed. If ix is missing, the global null hypothesis is used.}
  \item{tdp}{The threshold value for the intersection hypothesis which should be a proportion.}
}

\value{Returns the familywise error adjusted p-value.}

\author{Jelle Goeman.}

\references{
Meijer, R. J., Krebs, T. J., & Goeman, J. J. (2019). Hommel's procedure in linear time. Biometrical Journal, 61(1), 73-82.

Goeman, J., Meijer, R., Krebs, T., & Solari, A. (2016). Simultaneous control of all false discovery proportions in large-scale multiple hypothesis testing. arXiv preprint arXiv:1611.06739.
}

\examples{
  #Generate a vector of pvalues
  set.seed(102)
  m <- 10
  pvalues <- c(runif(0.5*m,0,0.02), runif(0.5*m,0,1))

  # First step: create a hommel object.
  hom <- hommel(pvalues, simes = TRUE)

  # Test any intersection hypothesis of interest against any chosen threshold
  localtest(hom, tdp=0.2)
  localtest(hom, 1:3, 0.8)

}
