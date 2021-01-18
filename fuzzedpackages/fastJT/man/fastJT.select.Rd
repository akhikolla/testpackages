\name{fastJT.select}
\alias{fastJT.select}
\docType{package}
\title{
  Conduct k-Fold Cross-Validation based on Jonckheere-Terpstra Test Statistics
}
\description{
  A method to conduct \code{k}-fold cross-validation based the Jonckheere-Terpstra test statistics for large numbers of dependent and independent variables, with optional multi-threaded execution. The data are divided into \code{k} subsets, one of which is withheld while the remaining \code{k}-1 subsets are used to test the features. The process is repeated \code{k} times, with each of the subsets withheld exactly once as the validation data. The calculation of the standardized test statistics employs the null variance equation as defined by Hollander and Wolfe (1999, eq. 6.19) to account for ties in the data.
}

\usage{
fastJT.select(Y, X, cvMesh=NULL, kFold=10L, selCrit=NULL, outTopN=15L, numThreads=1L)  	
}

\arguments{
  \item{Y}{A matrix of continuous values, representing dependent variables, e.g. marker levels or other observed values. Row names should be sample IDs, and column names should be variable names. Required.
}
  \item{X}{A matrix of integer values, representing independent variables, e.g. SNP counts or other classification features. Row names should be sample IDs, and column names should be feature IDs. Required.
}
  \item{cvMesh}{A user-defined function to specify how to separate the data into training and testing parts. The inputs of this function are a vector representing the sample IDs and \code{kFold}, an integer representing the number of folds for cross validation. The output of this function is a list of \code{kFold} vectors of sample IDs forming the testing subset for each fold. The default value is \code{NULL}, and if no function is specified, the data are partitioned sequentially into \code{kFold} equal sized subsets. Optional.
}
  \item{kFold}{An integer to indicate the number of folds. Optional. The default value is 10.
}
  \item{selCrit}{A user-defined function to specify the criteria for selecting the top features. The inputs of this function are \code{J}, the matrix of statistics resulting from \code{fastJT}, and \code{P}, the matrix of p-values from \code{pvalues(J)}. The output is a data frame containing the selected feature ID for each trait of interest. Optional. The default value is \code{NULL}, and if no function is specified, the features with the largest standardized Jonckheere-Terpstra test statistics are selected.}
  \item{outTopN}{An integer to indicate the number of top hits to be returned when \code{selCrit=NULL}. Unused if \code{selCrit!=NULL}. Optional. The default value is \code{15L}.
}
  \item{numThreads}{A integer to indicate the number of threads to be used in the computation. Optional. The default value is \code{1L} (sequential computation).
}
 
}

\note{
  Rows (samples) are assumed to be in the same order in \code{X} and \code{Y}.
}

\references{
    Hollander, M. and Wolfe, D. A. (1999) \emph{Nonparametric Statistical Methods}. New York: Wiley, 2nd edition.
}

\value{
  Three \code{lists} of length \code{kFold}
  \item{J}{A \code{list} of matrices of standardized Jonckheere-Terpstra test statistics, one for each cross validation.}
  \item{Pval}{A \code{list} of matrices of p-values, one for each cross validation.}
  \item{XIDs}{A \code{list} of matrices of the selected feature IDs, one for each cross validation.}
}

\seealso{
  \code{\link{fastJT}}
}

\examples{
# Generate dummy data	
num_sample  <- 100
num_marker  <- 10
num_feature <- 500
set.seed(12345)
Data <- matrix(rnorm(num_sample*num_marker), num_sample, num_marker)
Features <- matrix(rbinom(num_sample*num_feature, 2, 0.5), num_sample, num_feature)
colnames(Data) <- paste0("Var:",1:num_marker)
colnames(Features) <- paste0("Ftr:",1:num_feature)

res <- fastJT.select(Y=Data, X=Features, cvMesh=NULL, kFold=5, 
                     selCrit=NULL, outTopN=5, numThreads=1)
res 
}

