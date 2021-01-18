\name{quad}
\alias{quad}
\title{
~ Quadratic fitting function: quad ~
}
\description{
 Fits a set of paired points with a quadratic curve. Returns the quadratic set of points. Function only called by r_throwchart.
}
\usage{
quad(point,before_point,after_point,offSet)
}
\arguments{
 \item{point}{[numeric] or [integer]: A set of points between before and after points.}

 \item{before_point}{[numeric] or [integer]: The first point of the curve.}

 \item{after_point}{[numeric] or [integer]: The last point of the curve.} 

 \item{offSet}{[integer]: Single value offset for the graph.}
}
\value{
 returns the quadratic point equivalent.
}