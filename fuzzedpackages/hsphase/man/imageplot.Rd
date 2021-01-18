\name{imageplot}
\alias{imageplot}

\title{
Image Plot of Blocking Structure}
\description{
Create an imageplot of the blocking structure.}

\usage{
imageplot(x, title, rv = FALSE, ...)
}

\arguments{
 
\item{x}{\code{matrix} blocking structure (output of \link{bmh} or \link{hbp} functions)

}
 \item{title}{\code{character} title of imageplot

}
  \item{rv}{\code{logical}
reverse the colour
}
  \item{\dots}{

}
}

\details{
White indicates regions of unknown origin, red and blue correspond to the two sire strands.
}
\author{
This is a modified version of a function written by Chris Seidel. \cr
http://www.phaget4.org/R/image_matrix.html
}

\seealso{
\code{\link{bmh}} and \code{\link{aio}}
}
\examples{
genotype <- matrix(c(
0,2,1,1,1,
2,0,1,2,2,
2,2,1,0,2,
2,2,1,1,1,
0,0,2,1,0), ncol = 5, byrow = TRUE) # each row contains the SNP of individuals
imageplot(bmh(genotype))
}
\keyword{block}
\keyword{image}
