\docType{data}
\name{roster}
\alias{roster}

\title{Mutant counts from Rosche and Foster (2000)}

\description{
  Sample from Table 3 of Rosche and Foster (2000).
}

\usage{
  data(roster)
}

\format{
  A vector of positive integers with size 52.
}

\examples{
# are parameters significantly larger than 1 ?
flan.test(roster,alternative="greater",mutations0=1,fitness0=1)
}


\source{W. A. Rosche and P. L. Foster, Determining mutation rates in bacterial populations,  
\emph{Methods}, 20(1):4-17 (2000)}

\keyword{datasets}
