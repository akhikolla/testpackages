\name{Concentration set}

\alias{concentration}

\title{Concentration set calculation}

\description{Calculates the concentration set threshold. The concentration set is the subset of the p-values that contains all discoveries at confidence level 1-alpha.}

\usage{ 
  concentration (hommel, alpha=0.05) 
}

\arguments{
  \item{hommel}{The hommel object for the complete multiple testing problem.}
  \item{alpha}{The type I error allowed.}
}

\value{Returns a p-value. P-values larger than that value contain no discoveries in this data set at this level of alpha and may be disregarded.}

\author{Jelle Goeman.}

\references{
Goeman, J., Meijer, R., Krebs, T., & Solari, A. (2016). Simultaneous control of all false discovery proportions in large-scale multiple hypothesis testing. arXiv preprint arXiv:1611.06739.
}

\examples{
  #Generate a vector of pvalues
  set.seed(102)
  m <- 10
  pvalues <- c(runif(0.5*m,0,0.02), runif(0.5*m,0,1))
  
  # First step: create a hommel object.
  # Use Simes'test (simes = TRUE) or Hommel's robust test (simes = FALSE)
  hom <- hommel(pvalues, simes = TRUE)

  # Find the concentration set bound
  threshold <- concentration(hom)
  
  # Find the concentration set itself
  set <- pvalues <= threshold
  sum(set)
}