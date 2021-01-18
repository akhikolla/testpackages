\name{ohg}
\alias{ohg}
\title{
Matrix of Opposing Homozygotes 
}
\description{
Creates a matrix of opposing homozygotes from the genotype matrix.
}
\usage{
ohg(genotypeMatrix)
}
\arguments{
  \item{genotypeMatrix}{
{\code{matrix} genotype (Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data)
}
}
}


\value{
Returns a square matrix (sample X sample) with the pairwise counts of opposing homozygotes.
}

\note{
This function can be slow with a large data set. The fast version of this function will be available after publish of the related manuscript.
}
\author{
Ferdosi, M. H., & Boerner, V. (2014). A fast method for evaluating opposing homozygosity in large SNP data sets. Livestock Science.
}

\seealso{
 \code{\link{rpoh}}
}
\examples{
genotype <- matrix(c(        
  2,1,0,                  
  2,0,0,                  
  0,0,2                   
  ), byrow = TRUE, ncol = 3)    
  
ohg(genotype)             
}

\keyword{opposing-homozygote}
\keyword{snp}
\keyword{genoytpe}
