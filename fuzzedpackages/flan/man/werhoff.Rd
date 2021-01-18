\docType{data}
\name{werhoff}
\alias{werhoff}

\title{Rifampin-resistant mutant counts from Werngren and Hoffner (2003)}

\description{
  Samples from Table 1 p. 1522 of Werngren and Hoffner (2003). 
  
  The mean final number of cells is also given for each sample.
  
  The coefficient of variation of the final number of cells is the same for all samples.
}

\usage{
  data(werhoff)
}

\format{
  A list of 2 elements :
  \itemize{
    \item \code{cvfn} : the coefficient of variation of final numbers of cells, for all samples.
    \item \code{samples} : a list of 13 lists. Each list includes :
    \itemize{
      \item A sample of mutants counts named "W<index of the sample".
      \item \code{mfn} : the mean final number of cells.
    }
  }
}

\examples{
# coefficient of variation of final numbers
Wcvfn <- werhoff$cvfn
# tenth sample
W10 <- werhoff$samples[[10]]

# estimates without cvfn
without <- mutestim(W10$mc,mfn=W10$mfn)
# estimates with cvfn
with <- mutestim(W10$mc,mfn=W10$mfn,cvfn=Wcvfn)
rbind(without,with)
}

\source{J. Werngren and S. E. Hoffner: Drug susceptible Mycobacterium tuberculosis Beijing genotype does not develop mutation-conferred resistance to Rifampin at an elevated rate, 
\emph{J. Clin. Microbiol.}, 41(4):1520-1524 (2003)}
\keyword{datasets}
