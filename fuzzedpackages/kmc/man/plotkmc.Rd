\name{plotkmc2D}
\alias{plotkmc2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the contour plot of log-likelihood around the H0 ( dim=2).
}
\description{
Given a kmc object, this function will produce contour plot if there were two constraints.
}
\usage{
plotkmc2D(resultkmc,flist=list(f1=function(x){x},f2=function(x){x^2}),range0=c(0.2,3,20))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultkmc}{
      S3 Object of kmcS3.
%%     ~~Describe \code{x} here~~
}
\item{flist }{ 
    list of two functions,flist=list( f1=function( x) { x} ,f2=function( x) { x^2} )  
} 
\item{range0}{
	A vector that helps to determine the range of the contour plot, i.e (center[1]-range0[1], center[2]-range0[2]) to (center+range0[1], center[2]+range0[2]).
	The third parameter defines the number of grids would be used.
}
}

\value{
\item{X}{x.grid }
\item{Y}{y.grid }
\item{Z}{grid value }   
}

\author{
Yifan Yang
}

\examples{

x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) 
d <- c( 1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1) 

f<-function( x) { x-3.7} 

myfun5 <- function( x)  { 
 x^2-16.5
} 
# construnct g as a LIST!

g=list( f1=f,f2=myfun5) ; 
kmc.solve( x,d,g) ->re0;

#plotkmc2D(re0) ->ZZ; # run this to generate contour plot
#Advanced PLOT option using ggplot2: not run
#library(reshape2)
#volcano3d <- melt(ZZ$Z)
#names(volcano3d) <- c("x", "y", "z")

#volcano3d$x <- ZZ$X[volcano3d$x];
#volcano3d$y <- ZZ$Y[volcano3d$y];

#library(ggplot2)
#v <- ggplot(volcano3d, aes(x, y, z=z));
#v +geom_tile(aes(fill = z)) + stat_contour()+scale_fill_gradientn("Custom 
#Colours",colours=grey.colors(10));
#c("lightblue","blue","green","yellow","orange","red")
#X11();
#qplot(x, y, z = z, data = volcano3d, stat = "contour", geom = "path")
}


