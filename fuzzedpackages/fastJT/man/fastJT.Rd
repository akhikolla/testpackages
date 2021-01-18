\name{fastJT}
\alias{fastJT}
\docType{package}
\title{
  Compute the Jonckheere-Terpstra Test Statistics
}
\description{
  A method to compute the Jonckheere-Terpstra test statistics for large numbers of dependent and independent variables, with optional multi-threaded execution. The calculation of the standardized test statistics employs the null variance equation as defined by Hollander and Wolfe (1999, eq. 6.19) to account for ties in the data.
}

\usage{
fastJT(Y, X, outTopN=15L, numThreads=1L, standardized=TRUE)  	
}

\arguments{
  \item{Y}{A matrix of continuous values, representing dependent variables, e.g. marker levels or other observed values. Row names should be sample IDs, and column names should be variable names. Required.
}
  \item{X}{A matrix of integer values, representing independent variables, e.g. SNP counts or other classification features. Row names should be sample IDs, and column names should be feature IDs. Required.
}
  \item{outTopN}{An integer to indicate the number of top independent variables to be reported for each dependent variable, based on the standardized Jonckheere-Terpstra test statistics. Optional. The default value is \code{15L}. If set to \code{NA}, all results are returned.
}
  \item{numThreads}{A integer to indicate the number of threads used in the computation. Optional. The default value is \code{1L} (sequential computation).
}
  \item{standardized}{A boolean to specify whether to return standardized Jonckheere-Terpstra statistics (\code{TRUE}) or non-standardized statistics (\code{FALSE}). Optional. The default value is \code{TRUE}.
} 
}

\note{
  Rows (samples) are assumed to be in the same order in \code{X} and \code{Y}.
}

\references{
    Hollander, M. and Wolfe, D. A. (1999) \emph{Nonparametric Statistical Methods}. New York: Wiley, 2nd edition.
}

\value{
  A \code{list} with two objects
  \item{J}{A matrix of the standardized/non-standardized Jonckheere-Terpstra test statistics,
depending on the value of the \code{standardized} argument.
}
  \item{XIDs}{If \code{outTopN} was specified, this object is a matrix of the column IDs of \code{X} corresponding to the top standardized Jonckheere-Terpstra test statistics for each dependent variable. Otherwise this is a vector of the column IDs of \code{X}.
}
}


\examples{
# Generate dummy data	
num_sample <- 100
num_marker <- 10
num_SNP		 <- 500
set.seed(12345)
Mark <- matrix(rnorm(num_sample*num_marker), num_sample, num_marker)
Geno <- matrix(rbinom(num_sample*num_SNP, 2, 0.5), num_sample, num_SNP)
colnames(Mark) <- paste0("Mrk:",1:num_marker)
colnames(Geno) <- paste0("SNP:",1:num_SNP)

res <- fastJT(Y=Mark, X=Geno, outTopN=5)
res
res <- fastJT(Y=Mark, X=Geno, outTopN=NA)
head(res) 
}

