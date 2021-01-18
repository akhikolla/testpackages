\name{hss}
\alias{hss}

\title{
Half-sib Family Splitter}
\description{
Splits the dataset into half-sib family groups based on a pedigree.
}
\usage{
hss(pedigree, genotype, check = TRUE)
}

\arguments{
  \item{pedigree}{
 \code{matrix}  the pedigree matrix should contain at least two columns, the first column with the half-sib IDs
and the second column with the sires IDs
}
  \item{genotype}{
\code{matrix} genotype matrix with SNP ordered by mapping position in the columns. Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data
} 


\item{check}{
\code{logical} check the genotype file for the possible errors
}

}

\note{
Pedigree must have at least two columns with sample ids (Column 1) and sire ids (Column 2).
}
\details{
Only half-sib groups that have more than 3 individuals will be returned.
}
\value{
Returns a list of numeric matrices, each matrix is a half-sib family.
}
\examples{
# Please run demo(hsphase)
}



