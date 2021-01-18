\name{ohd}
\alias{ohd}
\docType{package}
\title{
Opposing Homozygote Detection}
\description{
Counts the number of opposing homozygotes for each animal that caused a heterozygus site in the sire.
}
\usage{
ohd(genotypeMatrix, unique_check = FALSE, SNPs = 6000)
}
\arguments{
  \item{genotypeMatrix}{
{\code{matrix} half-sib genotypes (one half-sib per row, with SNP ordered by mapping position in the columns. Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data)
}

}
\item{unique_check}{
{\code{logical} check if samples uniquely originate an opposing homozygote at a locus
}

}

 \item{SNPs}{
{\code{integer} number of SNP to use
}

}
}

\value{
Returns a vector with the number of heterozygous sites that each sample caused.
}
\note{
This function can be used to identify pedigree errors; i.e., the outliers.
}
\author{
This method is suggested by Bruce Tier <btier@une.edu.au> to identify pedigree errors.
}

\examples{
genotype <- matrix(c(        
  2,1,0,                  
  2,0,0,                  
  0,0,2                   
  ), byrow = TRUE, ncol = 3)    
  
ohd(genotype)             
}

\keyword{opposing homozygote}


