\docType{data}
\name{david}
\alias{david}

\title{Mutant counts and final numbers from H. L. David (1970)}

\description{
  Samples from Table 1 and 2 of H. L. David (1970).
  
Classes have been replaced by medians, 1000 for the last class.
}

\usage{
  data(david)
}

\format{
  A list of 11 lists named "D<index of the sample>". Each list includes :
  \itemize{
    \item \code{sample} A sample of mutants counts.
    \item \code{fn} : the mean final numbers of cells, or a sample of final numbers (only for the last one).
  }
  
  The i-th sample of the list includes the i-th column of Table 1, except the last sample which includes Table 2.
}

\examples{
# samples 1-10 only give mean final numbers
D1 <- david[[1]]; D1
# estimate mutation probability
mutestim(D1$mc,mfn=D1$mfn)

# sample 11 has individual final numbers
D11 <- david[[11]]; D11
# with mean final number
mutestim(D11$mc,mfn=mean(D11$fn))
# with true final numbers
mutestim(D11$mc,D11$fn)
}

\source{H. L. David, Probability distribution of drug-resistant mutants in unselected populations of 
\emph{Mycobacterium tuberculosis}. \emph{Appl. Microbiol.}, 20(5):810-814 (1970).
}
\keyword{datasets}
