\name{r_throwchart}
\alias{r_throwchart}
\title{
~ R graphics function: r_throwchart ~
}
\description{
 Used when interactive is false and creates a plot through R of this data visualisation method. 
}
\usage{
r_throwchart(before, after,xlim, ylim, col, lwd, offSet)
}
\arguments{
 \item{before}{[numeric] or [integer]: A (non-empty) vector of data values.}

 \item{after}{[numeric] or [integer]: A (non-empty) vector of data values.}

 \item{col}{[character]: A vector of hex code colours, by default "#123".}
 
 \item{lwd}{[integer]: Line width, a column of line widths, by default value is 2.5.}
 
 \item{xlim}{[numeric]: 2 value colum with x limits.}

 \item{ylim}{[numeric]: 2 value colum with y limits.}

 \item{offSet}{[integer]: Single value offset for the graph.}
}
\examples{
 if(interactive()){
    r_throwchart(tibble(c(1,2)),tibble(c(2,8)),c(0,0),c(0,0),col = "blue", c(1,5), offSet = 1)
 }
}