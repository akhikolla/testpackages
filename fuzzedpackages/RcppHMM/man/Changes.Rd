\name{Change Log}
\alias{Changes}
\title{Changes Made to Package RcppHMM}
\description{
This page contains a listing of recent changes made to the package.
}

\details{
\enumerate{
\item More examples were added to some functions. (November 2017)

\item Since there are different classes of HMMs and each of them with the same algorithms, a verification step was added to avoid memory leaks and variable compatibility. (May 2017)

\item The class of HMM with observations being modelled by a Gaussian Mixture Model (GHMM) was updated to have also a multivariate version (see \code{\link{initGHMM}}). (July 2017)

\item The emission matrix of the GHMM model was divided into two parameters: Mu and Sigma. Mu is now a 2D matrix with number of rows equal to the dimensionality of the observation vector and the number of columns equal to the number of hidden states. Sigma is now a 3D matrix with number of rows and columns equal to the the dimensionality of the observation vector and the number of slices equal to the number of hidden states (see \code{\link{initGHMM}}). (July 2017)

}}

\keyword{documentation}


