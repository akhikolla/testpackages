\name{True Discoveries}

\alias{discoveries}
\alias{fdp}
\alias{tdp}


\title{Confidence bounds for the number of discoveries}

\description{Calculates (1-alpha) lower confidence bounds for the number of false null hypotheses (discoveries) or for the true discovery proportion, and (1-alpha) upper confidence bounds for the false discovery proportions. These bounds are for arbitrary subsets, and are simultaneous over all subsets. They are valid even if the subset is chosen after seeing the data.}

\usage{ 
  discoveries (hommel, ix, incremental=FALSE, alpha=0.05) 
  fdp (hommel, ix, incremental=FALSE, alpha=0.05) 
  tdp (hommel, ix, incremental=FALSE, alpha=0.05) 
}

\arguments{
  \item{hommel}{The hommel object for the complete multiple testing problem.}
  \item{ix}{The selection of hypotheses considered. Any selector that works on the original vector of p-values (index, negative index, logical, names) is allowed. If missing, all hypotheses are selected.}
  \item{incremental}{If TRUE, the number of discoveries (of fdp or tdp) is reported for incremental sets formed by adding the indices from ix in order. If FALSE, reports the numbr of discoveries in the full set only.}
  \item{alpha}{The type I error allowed.}
}

\value{Returns a single integer (discoveries) or proportion (fdp, tdp), or a vector if incremental = TRUE. These represent the lower confidence bound (discoveries, tdp) or upper confidence bound (fdp). The confidence intervals are one-sided: the upper confidence bound for discoveries is the size of the set, and for tdp 1. The lower bound for fdp is 0.}

\author{Jelle Goeman, Rosa Meijer, Thijmen Krebs, Mitra Ebrahimpoor.}

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

  hom
  summary(hom)
  
  # Find lower confidence bounds 
  # for the number of false null hypotheses (discoveries)
  # In any subset of interest
  discoveries(hom, 1:2)
  
  # a bound for the false discovery proportion:
  set <- 1:5
  fdp(hom, set)
  
  # or incremental for sets ix = (1), (1,2), (1,2,3), (1,2,3,4) and (1,2,3,4,5)
  fdp(hom, set, incremental = TRUE)
  
  # use alpha=0.5 to obtain a conservative point estimate
  tdp(hom, set, alpha=0.5)
}