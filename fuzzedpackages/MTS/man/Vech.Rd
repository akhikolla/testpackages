\name{Vech}
\alias{Vech}
\title{Half-Stacking Vector of a Symmetric Matrix
}
\description{Obtain the half-stacking vector of a symmetric matrix
}
\usage{
Vech(mtx)
}
\arguments{
  \item{mtx}{A symmetric matrix
}
}
\details{Stacking a matrix into a vector using data on and below the diagonal.
}
\value{a vector consisting of stacked elements of a symmetric matrix
}
\author{Ruey S. Tsay
}
\examples{
m1=matrix(c(1:9),3,3)
m2=(m1+t(m1))/2
v1=Vech(m2)
}




