\name{FKS}
\alias{FKS}
\title{Fast Kernel Smoothing}
\usage{
FKS(dataf, Tout, kernel = c("uniform", "triangular", "Epanechnikov",
  "biweight", "triweight", "Gaussian"), m = 51, K = 20)
}
\arguments{
\item{dataf}{A set of functional data given by a \code{dataf} object that are to be smoothed.}

\item{Tout}{vector of values in the domain of the functions at which the
resulting smoothed function is evaluated}

\item{kernel}{Kernel used for smoothing. Admissible values are \code{uniform}, 
\code{triangular}, \code{Epanechnikov}, \code{biweight}, \code{triweight} and \code{Gaussian}.
By default, \code{uniform} is used.}

\item{m}{Number of points in the grid for choosing the cross-validated bandwidth.}

\item{K}{Performs \code{K}-fold cross-validation based on randomly shuffled data.}
}
\value{
A \code{dataf} object corresponding to \code{Tout} of smoothed functional values.
}
\description{
Produces a kernel smoothed version of a function based on
the vectors given in the input. Bandwidth is selected using cross-validation.
}
\details{
A vector of the same length as \code{Tout}
corresponding to the values of the
function produced using kernel smoothing, is provided. Bandwidth is selected using the
\code{K}-fold cross-validation of randomly shuffled input values.
}
\examples{
d = 10
T = sort(runif(d))
X = T^2+ rnorm(d,sd=.1)
Tout = seq(0,1,length=101)

plot(T,X)
dataf = list(list(args=T,vals=X))
data.sm = FKS(dataf,Tout,kernel="Epan")
lines(data.sm[[1]]$args,data.sm[[1]]$vals,col=2)

datafs = structure(list(dataf=dataf,labels=1:length(dataf)),class="functional")
plot(datafs)
points(T,X)
data.sms = structure(list(dataf=data.sm,labels=1:length(data.sm)),class="functional")
plot(data.sms)

n = 6
dataf = list()
for(i in 1:n) dataf[[i]] = list(args = T<-sort(runif(d)), vals = T^2 + rnorm(d,sd=.1))
data.sm = FKS(dataf,Tout,kernel="triweight")
data.sms = structure(list(dataf=data.sm,labels=1:length(data.sm)),class="functional")
plot(data.sms)
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{functional}
\keyword{kernel}
\keyword{smoothing}
