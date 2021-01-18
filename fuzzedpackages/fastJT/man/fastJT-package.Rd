\name{fastJT-package}
\alias{fastJT-package}
\alias{fastJT-package}
\docType{package}
\title{Efficient Jonckheere-Terpstra Test Statistics
}
\description{This \code{Rcpp}-based package implements highly efficient functions for the calculation of the Jonkheere-Terpstra statistic. It can be used for a variety of applications, including feature selection in machine learning problems, or to conduct genome-wide association studies (GWAS) with multiple quantitative phenotypes. The code leverages \code{OpenMP} directives for multi-core computing to reduce overall processing time.
}
\details{
\tabular{ll}{
Package: \tab fastJT\cr
Type: \tab Package\cr
Version: \tab 1.0.6\cr
Date: \tab 2020-11-04\cr
License: \tab GPL-3\cr
}
Please see the example function calls below, or refer to the individual function documentation or the included vignette for more information. The package vignette serves as a tutorial for using this package. The technical details are provided in the reference cited below. Specifically, the calculation of the standardized test statistics employs the null variance equation as defined by Hollander and Wolfe (1999, eq. 6.19) to account for ties in the data.
}

\author{
Jiaxing Lin, Alexander Sibley, Ivo Shterev, and Kouros Owzar

Maintainer: Alexander Sibley <dcibioinformatics@duke.edu>
}

\references{
    Hollander, M. and Wolfe, D. A. (1999) \emph{Nonparametric Statistical Methods}. New York: Wiley, 2nd edition.
}

\keyword{ package }
\seealso{
    \code{\link{Rcpp}}
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
   
res <- fastJT(Y=Data, X=Features, outTopN=15)
res
res <- fastJT.select(Y=Data, X=Features, cvMesh=NULL, kFold=5, 
                     selCrit=NULL, outTopN=5, numThreads=1)
res
}
