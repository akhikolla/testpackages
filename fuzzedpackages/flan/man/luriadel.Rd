\encoding{UTF-8}
\docType{data}
\name{luriadel}
\alias{luriadel}


\title{Resistant bacteria counts from Luria and Delbruck (1943)}

\description{
  Samples of mutants counts from Table 2 of Luria and Delbrück (1943).
 }

\usage{
  data(luriadel)
}

\format{
  A list of 3 samples of mutants counts, each named "LD<index of the sample>". 
  \itemize{
    \item The \code{LD1} sample includes the first four columns and the last one (experiments number 1, 10, 11, 15 and 21b) of table 2 p. 504.
    \item The \code{LD2} sample includes the fifth and the sixth columns (experiments number 16 and 17)
    \item The \code{LD3} sample includes the next to last column (experiment number 21a).
  }
  
}

\examples{
# 1st sample
mutestim(luriadel[[1]])

# compare 2nd and 3rd samples
flan.test(luriadel[2:3]) 
}


\source{D. E. Luria and M. Delbrück: Mutations  of bacteria from virus sensitivity to virus resistance, 
\emph{Genetics}, 28:491-511 (1943)
}
\keyword{datasets}
