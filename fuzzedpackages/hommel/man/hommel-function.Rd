\name{hommel}

\alias{hommel}

\title{Calculates adjusted p-values of Hommel's method efficiently.}

\description{Calculates adjusted p-values of Hommel's method efficiently.}

\usage{ hommel (p, simes = TRUE) }

\arguments{
  \item{p}{A vector of p-values.}
  \item{simes}{If \code{TRUE}, the Simes test is used (Hommel's method), if \code{FALSE} Hommel's robust test is used.}
}

\value{Returns a \code{\link{hommel-class}} object.}

\author{Jelle Goeman, Rosa Meijer, Thijmen Krebs.}

\references{
Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75(2), 383-386.

Meijer, R. J., Krebs, T. J., & Goeman, J. J. (2019). Hommel's procedure in linear time. Biometrical Journal, 61(1), 73-82.

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

  hom
  summary(hom)
  
  # retrieve familywise error adjusted p-values
  p.adjust(hom)
  
  # Find lower confidence bounds 
  # for the number of false null hypotheses (discoveries)
  # In any subset of interest
  discoveries(hom, 1:2)
}