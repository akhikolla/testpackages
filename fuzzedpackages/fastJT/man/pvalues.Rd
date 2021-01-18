\name{pvalues}
\alias{pvalues}
\docType{package}
\title{
	Compute P-values Based on Jonckheere-Terpstra Test Statistics
}
\description{
	Method to compute the p-values for results from \code{fastJT}. \code{fastJT} must be run with \code{standardized=TRUE} in order to use this function.
}

\usage{
pvalues(object)  	
}

\arguments{
  \item{object}{A \code{fastJT} object that is the return of method \code{fastJT}. Required.
}
}


\value{
  A matrix of p-values with the same dimensions as the standardized statistics from \code{fastJT}. 
}

\seealso{
  \code{\link{fastJT}}
}

\examples{
# Generate dummy data	
num_patient <- 100
num_marker 	<- 10
num_SNP     <- 500
set.seed(12345)
Mark <- matrix(rnorm(num_patient*num_marker), num_patient, num_marker)
Geno <- matrix(rbinom(num_patient*num_SNP, 2, 0.5), num_patient, num_SNP)
colnames(Mark) <- paste0("Mrk:",1:num_marker)
colnames(Geno) <- paste0("SNP:",1:num_SNP)

res <- fastJT(Y=Mark, X=Geno, outTopN=5)
pvalues(res)
res <- fastJT(Y=Mark, X=Geno, outTopN=NA)
pvalues(res)
}

