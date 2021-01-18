\name{getSigTests}
\alias{getSigTests}

\title{
 Extract Significant Test Results for a Given Alpha or Different Multiple Testing Corrections.
}

\description{
 This function finds the p-values less than or equal to \code{alpha}. Alpha is as given or corresponds to an adjusted alpha.
 Different multiple testing correction methods to adjust alpha are available.
}
\usage{
  getSigTests(pvalues, alpha=0.05, method="plain")
}

\arguments{
  \item{pvalues}{A matrix, vector or an object of class \code{gmw}.}
  \item{alpha}{Significance level.}
  \item{method}{Method of multiple testing adjustment.}
}

\details{
 The argument \code{pvalues} is the result of a \code{gmw} test run or a matrix/vector of p-values. If \code{pvalues} is in matrix shape (several test
 methods applied to the same data), then each row corresponds to a different test method and the columns to different variables.  The option \code{method}
 specifies the method of multiple testing correction. Typical options are \code{"plain"} for no correction, \code{"bonferroni"} for a Bonferroni correction,
 \code{"simes"} for an improved Bonferroni correction and \code{"BH"} for a Benjamini-Hochberg correction. Please note that the Simes method and the
 Benjamini-Hochberg correction lead to the same results.
 For permutation test results there is also the Westfall and Young method \code{"maxT"} available. In order to perform this correction the option \code{keepPM=TRUE}
 has to be set in the \code{gmw} call in order to keep the required permutation matrix.

 Additional there are the two options \code{"csD"} and \code{"csR"}. Those calculate for each alpha between 0 and \code{alpha} the difference (=\code{"csD"})
 or the ratio (=\code{"csR"}) between observed and expected rejections and report an optimal alpha (and corresponding test rejections) for which these
 criteria are maximal. Please keep in mind that this method might return unreasonable large 'optimal' cut-off points.

 See also the function \code{\link[gMWT]{rejectionPlot}} for more details. 
}
\value{
  A list object of class 're' with the values:
  (In case \code{pvalues} is a matrix the output is a list with length \code{nrow(X)} and each list item contains a single \code{re} object)
  \item{sigTests}{Position of the significant tests.}
  \item{sigPvalues}{P-values of the significant tests.}
  \item{pvalues}{The original \code{pvalues} object.}
  \item{method}{Chosen method.}
  \item{alpha}{Chosen alpha.}
  \item{multAlpha}{Adjusted / Optimal alpha.} 
  \item{inputN}{Rows of \code{pvalues}.}
}

\seealso{
\code{\link{rejectionPlot}}, \code{\link{gmw}}
}

\author{
  Daniel Fischer
}

\references{
  Benjamini, Y. and Hochberg, Y. (1995): Controlling the false discovery rate: a practical and powerful approach to multiple testing. J.R. Statist. Soc. B, 57(1):289 - 300.

  Simes, R. J. (1986): An improved bonferroni procedure for multiple tests of significance. Biometrika, 73:751 - 754.

  Westfall, P.H. and Young, S.S. (1993): Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment. Wiley, New York.
}

\examples{
  X <- matrix(c(rnorm(500,2,1),rnorm(600,2,1),rnorm(400,2.2,1)),byrow=TRUE, ncol=10)
  colnames(X) <- letters[1:10]
  g <- c(rep(1,50),rep(2,60),rep(3,40))
  test <- gmw(X,g,test="kw",type="external")
  cs1 <- getSigTests(test)
  cs1

  simData <- runif(1000,0,1)
  simData <- c(simData,runif(200,0,0.01))
  simData2 <- runif(1000,0,1)
  simData2 <- c(simData2,runif(200,0,0.01))

  simDataMat <- rbind(simData,simData2)

  getSigTests(simDataMat, method="bon")
  getSigTests(simData, method="bon")
  getSigTests(simData2, method="bon")
  getSigTests(simData, method="sim")
  getSigTests(pvalues=simData, method="bh",alpha=0.05)
  getSigTests(pvalues=simData, method="csD")
  getSigTests(pvalues=simData, method="csR")

  set.seed(731)
  X <- matrix(c(rnorm(50,2,1),rnorm(60,2,1),rnorm(40,2.2,1)),byrow=TRUE, ncol=10)
  colnames(X) <- letters[1:10]
  g <- c(rep(1,5),rep(2,6),rep(3,4))
  X[12:15,1] <- X[12:15,1] + 5

  # Keep the permutation matrix in order to perform the W&Y adjustment
    testPM <- gmw(X,g,test="kw",type="perm",keepPM=TRUE)
  # Apply the Westfall& Young adjustment
    getSigTests(testPM,method="maxT")

}

\keyword{methods}
 
