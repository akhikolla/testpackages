\name{throwchart}
\alias{throwchart}
\title{
~ Main function: throwchart ~
}
\description{
 Throwchart method is useful for visualising paired data, such as before/after data sets. Each pair of points are set on a horizontal axis and joined by a parabola. The height of the parabola is proportional to the difference: after- before = difference. If the difference is negative (after < before) then the curve is drawn under the axis.
}
\usage{
throwchart(before, after, col, id, lwd, xlim, ylim, offSet, webinteract)
}
\arguments{
 \item{before}{[numeric] or [integer]: A (non-empty) vector of data values.}

 \item{after}{[numeric] or [integer]: A (non-empty) vector of data values.}

 \item{col}{[character]: A vector of hex code colours, by default "#123".}

 \item{id}{[factor]: Column of string or number identifiers.}

 \item{lwd}{[integer]: Line width, a column of line widths, by default value is 2.5.}

 \item{xlim}{[numeric]: 2 value colum with x limits.}

 \item{ylim}{[numeric]: 2 value colum with y limits.}

 \item{offSet}{[integer]: Single value offset for the graph.}

 \item{webinteract}{[bool]: Is the function used in interactive mode?}
}
\examples{
  if(interactive()){
  	throwchart(c(1,2),c(2,8),c("#000","#F82"),id = c("id1","id2"),c(1,5))
  	throwchart(c(1,2),c(2,8), offSet = 1, webinteract=TRUE)
  }
  if(!interactive()){
        throwchart(c(1,2),c(2,8), offSet = 1, webinteract=FALSE)
        throwchart(c(1,2),c(2,1), webinteract=FALSE)
        throwchart(c(1,2),c(2,8),c("#000","#F00"),id = c("id1","id2"),c(1,5), webinteract=FALSE)
	n <- 10
	Avant <- rnorm(n)
	Apres <- Avant + rnorm(n) + 10
	throwchart(Avant, Apres, xlim =  c(-4,14), webinteract = FALSE)
	throwchart(Avant, Apres, offSet = 0, webinteract = FALSE)
	throwchart(Avant, Apres, offSet = 8, webinteract = FALSE)
  }
}
