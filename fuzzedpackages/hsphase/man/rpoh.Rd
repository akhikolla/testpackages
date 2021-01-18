\name{rpoh}
\alias{rpoh}
\title{
Reconstruct Pedigree Based on Matrix of Opposing Homozygotes}
\description{
Reconstructs a half-sib pedigree based on a matrix of opposing homozygotes.}
\usage{
rpoh(genotypeMatrix, oh, forwardVectorSize = 30, excludeFP = TRUE, nsap = 3,
maxRec = 15, intercept = 26.3415, coefficient = 77.3171, snpnooh, method, maxsnpnooh)
}
\arguments{
  \item{genotypeMatrix}{
\code{matrix} genotype of one chromosome (data should be numeric. Use 0, 1 and 2 for respectively AA, AB and BB. Use 9 for missing data)
}
\item{oh}{
  \code{integer} Opposing homozygotes matrix (Output of \code{\link{ohg}}) }


  \item{forwardVectorSize}{
   \code{integer} number of heterozygous sites used to validate recombination events or check for genotyping errors}


\item{excludeFP}{
  \code{logical} excludes SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors}


\item{nsap}{
  \code{integer} number of SNP per block to validate recombinations}

\item{maxRec}{
  \code{integer} maximum number of expected recombinations per individual}
  
\item{intercept}{
  \code{integer} intercept of fitted model}
  
\item{coefficient}{
  \code{integer} coefficient of fitted model}
  
  \item{snpnooh}{
  \code{integer} number of SNPs used to create \emph{oh} matrix (this number must be divided by 1000)}
  
\item{method}{
  \code{character} pedigree reconstruction method}
  
\item{maxsnpnooh}{
  \code{numeric} the maximum number of allowing opposing homozygote in a half-sib family
}
}

\details{
Four methods \emph{simple}, \emph{recombinations}, \emph{calus} and \emph{manual} can be
utilized to reconstruct the pedigree. 
\cr\cr
The following examples show the arguments require for each method.

pedigree1 <- rpoh(oh = oh, snpnooh = 732, method = "simple") \cr
pedigree2 <- rpoh(genotypeMatrix = genotypeChr1, oh = ohg(genotype), maxRec = 10 , method = "recombinations")  \cr
pedigree3 <- rpoh(genotypeMatrix = genotype, oh = oh, method = "calus") \cr
pedigree4 <- rpoh(oh = oh, maxsnpnooh = 31662, method = "manual")
 \cr
}

\value{
Returns a data frame with two columns, the first column is animals' ID and the second column is sire identifiers (randomly generated).
}
\note{
Method can be \emph{recombinations}, \emph{simple}, \emph{calus} or \emph{manual}. Please refer to vignette for more information.\cr

The sire genotype should be removed before using this function utilizing  \code{\link{pogc}} function.
}



\seealso{
 \code{\link{bmh}} and \code{\link{recombinations}} 
}
\examples{
# Please run demo(hsphase)
}

\keyword{snp}
\keyword{opposing-homozygote}
\keyword{pedigree}
\keyword{reconstruction}