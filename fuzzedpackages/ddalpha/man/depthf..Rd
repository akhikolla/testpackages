\name{depthf.}
\alias{depthf.}
\title{
Calculate Functional Depth
}
\description{
Calculates the depth of functions w.r.t. a functional data set.

The detailed descriptions are found in the corresponding topics.
}
\usage{
depthf.(datafA, datafB, notion, ...)

## Adjusted band depth
# depthf.ABD(datafA, datafB, range = NULL, d = 101, norm = c("C", "L2"), 
# J = 2, K = 1)

## Band depth
# depthf.BD(datafA, datafB, range = NULL, d = 101)

## Univariate integrated and infimal depth
# depthf.fd1(datafA, datafB, range = NULL, d = 101, order = 1, approx = 0)

## Bivariate integrated and infimal depth
# depthf.fd2(datafA, datafB, range = NULL, d = 101)

## h-mode depth
# depthf.hM(datafA, datafB, range = NULL, d = 101, norm = c("C", "L2"),
#  q = 0.2)

## Bivariate h-mode depth
# depthf.hM2(datafA, datafB, range = NULL, d = 101, q = 0.2)

## Half-region depth
# depthf.HR(datafA, datafB, range = NULL, d = 101)

## Univariate random projection depths
# depthf.RP1(datafA, datafB, range = NULL, d = 101, nproj = 50, nproj2 = 5)

# Bivariate random projection depths
# depthf.RP2(datafA, datafB, range = NULL, d = 101, nproj = 51)
}
\arguments{
  \item{datafA}{
Functions whose depth is computed, represented by a \code{dataf} object of their arguments
and functional values.
}
  \item{datafB}{
Random sample functions with respect to which the depth of \code{datafA} is computed. 
\code{datafB} is represented by a \code{dataf} object of their arguments
 and functional values.
}
  \item{notion}{
The name of the depth notion (shall also work with a user-defined depth function named \code{"depthf.<name>"}).
}
  \item{\dots}{
Additional parameters passed to the depth functions.
}
}

\seealso{

\code{\link{depthf.ABD}}

\code{\link{depthf.BD}}

\code{\link{depthf.fd1}}

\code{\link{depthf.fd2}}

\code{\link{depthf.hM}}

\code{\link{depthf.hM2}}

\code{\link{depthf.HR}}

\code{\link{depthf.RP1}}

\code{\link{depthf.RP2}}

}
\value{
Numerical vector of depths, one for each function in \code{datafA}; or one depth value if \code{datafA} is a single function.
}
\examples{
# real data example
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]
                
depthf.(datafA, datafB, notion = "HR")

dataf2A = derivatives.est(datafA,deriv=c(0,1))
dataf2B = derivatives.est(datafB,deriv=c(0,1))

depthf.(dataf2A, dataf2B, notion = "fd2")
}
\keyword{ robust }
\keyword{ functional }
\keyword{ nonparametric }
