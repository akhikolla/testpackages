\name{recombinations}
\alias{recombinations}
\title{
Recombination Number
}
\description{
Counts the number of recombinations for each individual.
}
\usage{
recombinations(blockMatrix)
}
\arguments{
  \item{blockMatrix}{
\code{matrix} block structure (Output of \code{\link{bmh}})
}
}

\value{
Returns a vector of recombinations. The number of elements in this vector is equal to the number of individuals,  
i.e. each element holds the number of recombinations identified for each sample.
}



\seealso{
 \code{\link{bmh}
}}
\examples{
genotype <- matrix(c(          
  2,1,0,0,                    
  2,0,2,2,                    
  0,0,2,2,
  0,2,0,0                     
  ), byrow = TRUE, ncol = 4)      
  
recombinations(bmh(genotype))             
}

\keyword{recombination}
\keyword{pedigree}
\keyword{error}
\keyword{block}
