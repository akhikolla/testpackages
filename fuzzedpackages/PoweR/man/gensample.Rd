\name{gensample}
\alias{gensample}
\alias{rlaw}
\title{Generate random samples from a law added in the package.}
\description{Generate random samples from a law added in the package as a C function.}
\usage{
gensample(law.index,n,law.pars = NULL,check = TRUE, center=FALSE, scale=FALSE)
}
\arguments{
  \item{law.index }{law index as given by function
    \code{\link{getindex}}.}
  \item{n }{number of observations to generate.}
  \item{law.pars }{vector of parameters for the law. The length of this
    parameter should not exceed 4.}
  \item{check }{logical. If \code{TRUE}, we check if \code{law.index} belongs to 
    the list of laws. If \code{FALSE}, we pass on this verification, this will
	reduce the simulation time.}
  \item{center}{Logical. Should we center the data generated}
  \item{scale}{Logical. Should we center the data generated}
}
\value{A list containing the random sample and the vector of parameters used for the chosen law.}
\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}

\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{checklaw}}}
\examples{
# This is good to check if the generator of the given law has been well coded.

res <- gensample(2,10000,law.pars=c(-5,2),check=TRUE)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)

# See function checklaw() in this package.
hist(gensample(2,10000,law.pars=c(0,1),check=TRUE)$sample,prob=TRUE,breaks=100,main="Density
histogram of the N(0,1) distribution")
curve(dnorm(x),add=TRUE,col="blue")

}
\keyword{distribution}
\keyword{datagen}
