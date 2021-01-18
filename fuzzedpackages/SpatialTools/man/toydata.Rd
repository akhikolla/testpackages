\name{toydata} 

\alias{toydata}

\title{
	A toy data set for use in examples.
} 

\description{
	A list containing \code{X}, a 50 x 3 design matrix, \code{y}, a vector of length 50 of observed responses, \code{V}, a 50 x 50 covariance matrix for the observed data, \code{Xp}, a 121 x 3 design matrix for the predicted responses, \code{Vp}, the 121 x 121 covariance matrix of the predicted responses, \code{Vop}, the 50 x 121 covariance matrix between the observed responses and the predicted responses, \code{coords}, a 50 x 2 matrix containing the sites of the 50 observed responses, and \code{pcoords}, a 121 x 2 matrix containing the 121 sites for the predicted responses (a 11 x 11 regular grid over the domain [0, 1]x[0, 1]).
} 

\usage{
	data(toydata)
} 

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\examples{
	data(toydata)
}