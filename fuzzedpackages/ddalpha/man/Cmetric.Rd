\name{Cmetric}
\alias{Cmetric}
\title{Fast Computation of the Uniform Metric for Sets of Functional Data}
\usage{
Cmetric(A, B)
}
\arguments{
\item{A}{Functions of the first set, represented by a matrix of their functional values of 
size \code{m*d}. \code{m} stands for the number of functions, \code{d}
is the number of the equi-distant points in the domain of the data at which the functional
values of the \code{m} functions are evaluated.}

\item{B}{Functions of the second set, represented by a matrix of their functional values of 
size \code{n*d}. \code{n} stands for the number of functions, \code{d}
is the number of the equi-distant points in the domain of the data at which the functional
values of the \code{n} functions are evaluated. The grid of observation points for the 
functions \code{A} and \code{B} must be the same.}
}
\value{
A symmetric matrix of the distances of the functions of size \code{m*n}.
}
\description{
Returns the matrix of \eqn{C} (uniform) distances between two sets of functional data.
}
\details{
For two sets of functional data of sizes \code{m} and \code{n}
represented by matrices of their functional values, 
this function returns the symmetric matrix of size \code{m*n} whose entry in the
\code{i}-th row and \code{j}-th column is the approximated \eqn{C} (uniform) distance of the 
\code{i}-th function from the first set, and the \code{j}-th function from the second set.
This function is utilized in the computation of the h-mode depth.
}
\examples{
datapop = dataf2rawfd(dataf.population()$dataf,range=c(1950,2015),d=66)
A = datapop[1:20,]
B = datapop[21:50,]
Cmetric(A,B)

}
\seealso{
\code{\link{depthf.hM}}

\code{\link{dataf2rawfd}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{functional}
\keyword{metric}
