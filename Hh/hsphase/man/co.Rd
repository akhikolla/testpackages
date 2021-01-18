% 
% TODO: Write Documentation 
% 

\name{co}
\alias{co}
\docType{package}
\title{
Crossover Detection}
\description{
Detect all possible crossover events.}
\usage{
co(genotypeMatrix)
}
\arguments{
  \item{genotypeMatrix}{
{\code{matrix} half-sib genotypes  (one half-sib per row, with SNP ordered by mapping position in the columns. Data should be numeric. Use 0, 1 and 2 for respectively AA, AB and BB. Use 9 for missing data)
}

}
}

\value{
Returns a matrix with the number of crossover events for each site.
}


\examples{
genotype <- matrix(c(           # Define a Half-sib Genotype Matrix
  2,1,0,                        # Individual 1
  2,0,2,                        # Individual 2
  0,0,2                         # Individual 3
  ), byrow = TRUE, ncol = 3)    # There are 3 individuals with three SNPs
  
co(genotype)               

}

\keyword{crossover}




