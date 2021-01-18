\name{ISOpure.util.repmat}
\alias{ISOpure.util.repmat}
\title{Tiles matrix horizontally or vertically}
\description{Tiles matrix horizontally or vertically in the same way as the Matlab repmat command}
\usage{
ISOpure.util.repmat(a, n, m)
}
\arguments{
  \item{a}{A matrix}
  \item{n}{Number of times the matrix should be tiled horizontally}
  \item{m}{number of times the matrix should be tiled vertically}
}
\value{A matrix which has replicated and tiled the input matrix a by n rows and m columns}
\author{Catalina Anghel, Ohloh (now Black Duck Open Hub)}
\examples{
x <- matrix(runif(6), 3, 2)
x
#           [,1]      [,2]
# [1,] 0.5167029 0.7543404
# [2,] 0.9064936 0.4316977
# [3,] 0.3256870 0.5310625
ISOpure.util.repmat(x, 1, 2)
#           [,1]      [,2]      [,3]      [,4]
# [1,] 0.5167029 0.7543404 0.5167029 0.7543404
# [2,] 0.9064936 0.4316977 0.9064936 0.4316977
# [3,] 0.3256870 0.5310625 0.3256870 0.5310625
ISOpure.util.repmat(x, 2, 1)
#           [,1]      [,2]
# [1,] 0.5167029 0.7543404
# [2,] 0.9064936 0.4316977
# [3,] 0.3256870 0.5310625
# [4,] 0.5167029 0.7543404
# [5,] 0.9064936 0.4316977
# [6,] 0.3256870 0.5310625
ISOpure.util.repmat(x, 2, 3)
#           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
# [1,] 0.5167029 0.7543404 0.5167029 0.7543404 0.5167029 0.7543404
# [2,] 0.9064936 0.4316977 0.9064936 0.4316977 0.9064936 0.4316977
# [3,] 0.3256870 0.5310625 0.3256870 0.5310625 0.3256870 0.5310625
# [4,] 0.5167029 0.7543404 0.5167029 0.7543404 0.5167029 0.7543404
# [5,] 0.9064936 0.4316977 0.9064936 0.4316977 0.9064936 0.4316977
# [6,] 0.3256870 0.5310625 0.3256870 0.5310625 0.3256870 0.5310625
}
\keyword{manip}
\keyword{array}