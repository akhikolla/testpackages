\name{get.contours} 

\alias{get.contours}

\title{
	Extracts coordinates from contourLines function
} 

\description{
	Takes contours of \code{contourLines} function and extracts the associated coordinates.
} 

\usage{
	get.contours(x)
} 

\arguments{
	\item{x}{A list returned by the \code{contourLines} function.} 
} 

\value{
Returns a 2-column matrix containing the coordinates making up the contours in contours.list.
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{contour}
\keyword{contour lines}
\keyword{level curve}


\seealso{contourLines, contour}

\examples{
data(volcano)
x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
cL <- contourLines(x, y, volcano)
out <- get.contours(cL)
}