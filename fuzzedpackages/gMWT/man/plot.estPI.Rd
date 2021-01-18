\name{plot.estPI}
\docType{methods}
\alias{plot.estPI}
\alias{plot-method}
\alias{plot,estPI-method}

\title{Plot an estPI Object}
\description{
  Plots an \code{estPI} object.
}
\usage{
 \method{plot}{estPI}(x,col="black",highlight,hlCol="red",pch=20,zoom=FALSE,...)
}

\arguments{
  \item{x}{Object of class \code{estPI}.}
  \item{col}{Vector of colors of the scatterplot.}
  \item{highlight}{Vector with positions, which shall be marked in special color.}
  \item{hlCol}{Color of highlighted spots.}
  \item{pch}{Dot type of the plot.}
  \item{zoom}{Logical, shall the plots be zoomed to interesting areas?}
  \item{...}{Additional parameters.}
}

\details{
This function plots the probabilistic indices from an \code{estPI} object. Additional graphical
parameters are \code{col}, that specifies the colors of the points used in the plot. Furthermore there
is the option to highlight the probabilistic indices for specified variables. Those are defined with
the \code{highlight} option and the \code{hlCol} color. One possible application of this is to highlight
only the significant test results, or e.g. in the genetic context only those genes, which are located
on a certain chromosome.

The boolean option \code{zoom} is meant for plotting the graph in a different scale and it adjusts the
limits of the x and y axes to the minimum and maximum of the corresponding probabilistic index. 
}

\examples{
  X <- c(sample(15))
  g <- c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,4)
  res <- estPI(X,g,type="single")
  plot(res)

  X <- matrix(c(rnorm(5000,1.5,2),rnorm(6000,2,2),rnorm(4000,3.5,1),rnorm(3000,2.5,1)),
           byrow=TRUE, ncol=100)
  colnames(X) <- letters[1:100]
  g <- c(rep(1,50),rep(2,60),rep(3,40),rep(4,30))

  res <- estPI(X,g,type="single")
  plot(res)
}



\author{Daniel Fischer}

\keyword{methods}
\keyword{plot}
