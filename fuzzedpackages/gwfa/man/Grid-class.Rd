\docType{class}
\name{Grid-class}
\alias{Grid-class}
\title{Grid class}
\alias{[,Grid-method}
\alias{[<-,Grid-method}

\description{
A S4 class to define a grid, the output class of gwfa function.
}

\value{
an object of class "Grid"
}

\section{Slots}{
\describe{
\item{cell_size}{Cell size of the grid}
\item{bandwith}{Radius of the Kernel Density Estimator}
\item{radius}{Radius of the multiscale analysis}
\item{q}{The vector of order of the Renyi entropy}
\item{sample_size}{size of the sample for calculate the average of the masses for the different radius and Renyi entropy orders \eqn{M^q(R)}
}
}}

\section{Methods}{
\describe{
\item{\code{[}}{Gets the value of an object}
\item{\code{[<-}}{Sets the value of an object}
}
}