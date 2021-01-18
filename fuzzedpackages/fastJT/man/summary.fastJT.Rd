\name{summary.fastJT}
\alias{summary.fastJT}
\docType{package}
\title{
  Summarize Jonckheere-Terpstra Test Statistics and P-Values
}
\description{
  Summary method for fastJT results.
}

\usage{
\method{summary}{fastJT}(object, Y2Print=1:10, X2Print=1:10, printP=TRUE, 
               outTopN=NA, subObj=FALSE, \dots)	
}

\arguments{
  \item{object}{A \code{fastJT} object that is the return of method \code{fastJT}. Required.
}
  \item{Y2Print}{Either a numeric or character vector that indicates the desired dependent variables to print. The default is 1:10. Set to \code{NA} to print full results. Optional.
}
  \item{X2Print}{Either a numeric or character vector that indicates the desired independent variables to print. If \code{outTopN=NA} in function \code{fastJT}, i.e., the results are not sorted, both numeric and character vectors can be used to set the print range. The default range is 1:10. Set to \code{NA} to print full results. If \code{outTopN!=NA} in function \code{fastJT}, the range of \code{X2Print} refers to the range of top normalized statistics computed in \code{fastJT}, and only numeric vectors can be used. Optional.
}
  \item{printP}{A boolean indicating whether to print the p-values (\code{TRUE}) or the standardized statistics (\code{FALSE}). The default value is \code{TRUE}. Optional.
}
  \item{outTopN}{An integer specifying the number of top hits to print in the summary, if the statistics were not sorted during the evaluation of the statistics (i.e., if \code{outTopN=NA} in function \code{fastJT}). Optional. The default value is \code{NA}.
}
  \item{subObj}{A boolean indicating whether to return a \code{fastJT} object subset per the requested summary (\code{TRUE}). Optional. The default value is \code{FALSE} (nothing is returned).
}
  \item{...}{Additional arguments affecting the summary produced.
}
}

\value{
  If \code{subObj=TRUE}, this method returns a \code{fastJT} object matching the subset of statistics from \code{object} that are being printed. For example, if \code{object} is not sorted by top hits, \code{summary(object, outTopN=10, subObj=TRUE)} will print the summary and return a subset of \code{object} that contains only the top 10 independent variables for each dependent variable. If \code{subObj=FALSE}, the summary is printed but no values are returned.
}
\note{
  This function prints a matrix or paired columns of independent variable IDs and statistics/p-values to the log.
}
\seealso{
  \code{\link{fastJT}}, \code{\link{pvalues}}
}

\examples{
# Generate dummy data	
num_patient <- 100
num_marker 	<- 10
num_SNP     <- 500
set.seed(12345);
Mark <- matrix(rnorm(num_patient*num_marker), num_patient, num_marker)
Geno <- matrix(rbinom(num_patient*num_SNP, 2, 0.5), num_patient, num_SNP)
colnames(Mark) <- paste0("Mrk:",1:num_marker)
colnames(Geno) <- paste0("SNP:",1:num_SNP)

res <- fastJT(Y=Mark, X=Geno, outTopN=5)
summary(res, Y2Print=c("Mrk:1","Mrk:2"), X2Print=1:5, printP=FALSE)
summary(res, Y2Print=NA, X2Print=1:5, printP=FALSE)

res <- fastJT(Y=Mark, X=Geno, outTopN=NA)
summary(res, Y2Print=1:10, X2Print=1:10, printP=TRUE)
summary(res, Y2Print=c("Mrk:1","Mrk:2"), X2Print=c("SNP:1","SNP:2"), printP=TRUE) 

res <- fastJT(Y=Mark, X=Geno, outTopN=NA, standardized=FALSE)
summary(res, outTopN=10)
}

