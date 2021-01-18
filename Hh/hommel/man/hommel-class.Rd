\name{hommel class}
\docType{class}
\alias{hommel-class}
\alias{p.adjust,hommel-method}
\alias{summary,hommel-method}
\alias{show,hommel-method}

\title{Class "hommel" for storing the result of Hommel's multiple testing method.}

\description{The class hommel is the output of a call to \code{\link{hommel}}. It stores the adjusted p-values of all elementary hypotheses, as well as information needed to quickly calculate other results about the underlying closed testing procedure.}

\section{Slots}{
  These slots are not meant to be directly accessed by the user.
  \describe{
    \item{\code{p}:}{Object of class "numeric". Stores original (unadjusted) p-values for all hypotheses.}
    \item{\code{sorter}:}{Object of class "integer". Stores the permutation to sort p.}
    \item{\code{jumpalpha}:}{Object of class "numeric". Stores the alpha-values where the function h jumps.}
    \item{\code{adjusted}:}{Object of class "numeric". Stores adjusted p-values for all elementary hypotheses.}
    \item{\code{simesfactor}:}{Object of class "numeric". The denominator of the local test that is used.}
    \item{\code{simes}:}{Object of class "logical". Is TRUE if the adjusted p-values are calculated based on Simes' test that assumes the Simes inequality. If FALSE if they are based on Hommel's robust test.}
  }
}


\section{Methods}{
  \describe{
    \item{p.adjust}{(hommel): Retrieves adjusted pvalues for all elementary hypotheses.}
    \item{summary}{(hommel, alpha=0.05): Summarizes the results of hommel's method for the specified level of alpha.}
  }
}

\seealso{\link{select}, \link{concentration}.}

\author{Jelle Goeman, Rosa Meijer, Thijmen Krebs.}

\references{
Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75(2), 383-386.

Meijer, R. J., Krebs, T. J., & Goeman, J. J. (2019). Hommel's procedure in linear time. Biometrical Journal, 61(1), 73-82.

Goeman, J., Meijer, R., Krebs, T., & Solari, A. (2016). Simultaneous control of all false discovery proportions in large-scale multiple hypothesis testing. arXiv preprint arXiv:1611.06739.
}

\keyword{methods}

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