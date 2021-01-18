\name{tt.nodes}
\alias{tt.nodes}
\docType{data}
\title{
   Host list with there status over the time.
}
\description{
   a data.frame of all the hosts identify by there ID.
  Three colums is use ID, status and time
}
\usage{
data("transmissiontree")
}
\format{
  A data frame with 47 observations on the following 3 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{status}}{a character vector}
    \item{\code{time}}{a character vector}
  }
}
\examples{
data(transmissiontree)
print(tt.nodes)
}
\keyword{datasets}
