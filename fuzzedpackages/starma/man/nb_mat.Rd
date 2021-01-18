\name{nb_mat}
\docType{data}
\alias{nb_mat}
\alias{dlist}
\alias{klist}
\alias{blist}
\title{ Neighbourhood weight matrices for France's 94 departments }
\description{
This data file provides three neighbourhoods for the 94 metropolitan French departments:

- \code{dlist}: distance-based neighbourhoods; two departments are considered neighbours if their centroids are within range of 100km.

- \code{klist}: four closest neighbours; each department is connected to its four closest neighbours, the distance being calculated between centroids.

- \code{blist}: common border neighbours; two departments are considered neighbours if they share a border.

These neighbourhoods are designed to be used within the \code{\link{starma-package}}.
First element is the identity matrix (0-th order neighbours).
Second element is the common border contingency matrix of the department (1-st order neighbours).
Elements three to five are the weight matrices lagged from the previous one (2-nd to 4-th order neighbours).

They have been computed used the package \code{spdep} and its functions \code{readShapePoly}, \code{poly2nb} and \code{nblag}.
}
\usage{nb_mat}
\format{
Three lists of 5 weight matrices, of dimension 94x94
}
\keyword{datasets}
\keyword{neighbourhood}
\keyword{weight matrix}
\keyword{dlist}
\keyword{klist}
\keyword{blist}
