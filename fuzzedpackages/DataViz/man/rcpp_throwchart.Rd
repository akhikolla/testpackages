\name{rcpp_throwchart}
\alias{rcpp_throwchart}
\title{
  ~ C++ called function: rcpp_throwchart ~
}
\description{
  This function takes the inputs from \link{throwchart}, and writes the data in Json array, then this function calls a windows cmd function to open an index.html in the extdata.
}
\usage{
rcpp_throwchart(before, after ,col, id, lwd, xlim, ylim, offSet, path)	
}
\arguments{
  \item{before}{A number column}

  \item{after}{A number column}
  
  \item{col}{A hex code colour colum has to be format "#000"}
  
  \item{id}{An id has to be string}
  
  \item{lwd}{A number for the line width best between 1-5)}
  
  \item{xlim}{[numeric]: 2 value colum with x limits.}
  
  \item{ylim}{[numeric]: 2 value colum with y limits.}
  
  \item{offSet}{[integer]: Single value for the graph offset.}

  \item{path}{The path for the library}
}
\examples{
 if(interactive()){
 before = tibble(c(1,2))
 after = tibble(c(2,8))
 col = tibble(c("#000","#F82"))
 id = tibble(c("",""))
 lwd = tibble(c(1,5))
 xlim = tibble(c(0,0))
 rcpp_throwchart(before,after,col,id,lwd,xlim,0,path.package("DataViz"))}
}
