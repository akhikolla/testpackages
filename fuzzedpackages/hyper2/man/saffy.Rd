\name{saffy}
\alias{saffy}
\title{Choice matrix functionality}
\description{
Converts a matrix of restricted choices into a likelihood function.  It
is named for Saffron O'Neill.
}
\usage{
saffy(M)
}

\arguments{
  \item{M}{Matrix}
}
\details{
The canonical example is Saffron's climate change dataset, documented at
\code{icons}.  Function \code{saffy()} returns the appropriate
likelihood function for the dataset.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{icons}}}
\examples{

data("oneill")
icons == saffy(icons_matrix)  # should be TRUE

}