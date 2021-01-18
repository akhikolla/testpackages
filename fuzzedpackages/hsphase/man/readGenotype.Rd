% 
% TODO: Write Documentation 
% 

\name{readGenotype}
\alias{readGenotype}

\title{
Read and Check the Genotype File}
\description{
This function reads and checks genotype files.
}
\usage{
readGenotype(genotypePath, separatorGenotype = " ", check = TRUE)
}

\arguments{

  \item{genotypePath}{
\code{character} genotype path (animals (rows) and SNP (columns), SNP should be coded as 0, 1 and 2 for respectively AA, AB and BB.
Use 9 for missing data. please refer to vignette for more information)
}

  \item{separatorGenotype}{
 \code{character} separator character for genotype
}

\item{check}{
\code{logical} check the genotype file for possible errors
}

}

\note{
Please refer to vignette for more information.
}


\value{
Returns the genotype matrix.
}
\examples{
# A comprehensive demo and example dataset is available from
# http://www-personal.une.edu.au/~cgondro2/hsphase.html
}
