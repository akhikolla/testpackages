\name{DataViz-package}
\alias{DataViz-package}
\alias{DataViz}
\docType{package}

\title{
 ~ Overview: package DataViz ~
}

\description{Various data vizualisation methods.}

\details{
 Data Visualisation is the art of graphically representing data.
 There are numerous data visualisation methods, but they aren't always relevant -and sometimes less informative than basic representations-.
 Moreover they are often created by programmers in various computer languages and the code being seldom available.

 DataViz is a package aiming to give access to Data Visualisation methods that are relevant from the statistician's point of view.

 The 3 first methods to be implemented are \link{throwchart}, Gravity Bubble Chart (V0.3, june 2019) and XXX (V0.4, septembre 2019)
}

\author{
Timothy Bell \email{horia.yeb@gmail.com}
Christophe Genolini \email{christophe.genolini@u-paris10.fr}

Maintainer: Timothy Bell <horia.yeb@gmail.com>
}

\references{
 Inspired from http://tiffanyfrance.com/data-is-beautiful/19-01/
}

\seealso{
  \link{throwchart}
  \link{forcelayout}
}

\examples{
 if(interactive()){
   throwchart(c(1,2),c(2,8),c("#000","#F82"),id = c("id1","id2"),c(1,5))
   throwchart(c(1,2),c(2,8))
}
if(!interactive()){
   throwchart(c(1,2),c(2,8), offSet = 1, webinteract=FALSE)
   throwchart(c(1,2),c(2,1), webinteract=FALSE)
   throwchart(c(1,2),c(2,8),c("#000","#F00"),c(1,5), webinteract=FALSE)
}
}
\keyword{ DataViz }
