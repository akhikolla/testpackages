\name{KernSmoothIRT-package}
\alias{KernSmoothIRT}
\docType{package}
\title{KernSmoothIRT Package}
\description{
Graphical analysis of multiple choice test and questionnaire data within a nonparametric approach. Fits item and option characteristic curves using kernel smoothing techniques and provides a variety of exploratory plots.
}
\details{
\tabular{ll}{
Package: \tab KernSmoothIRT\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-29\cr
Requires: \tab plotrix \cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

The function \code{ksIRT} creates a \code{ksIRT} object using kernel smoothing. 
Plotting the \code{ksIRT} object creates a variety of exploratory plots, including: item characteristic curves, option characteristic curves, density plots, expected value plots, and a probability simplex plot for the top 3 highest probability options of each item. 

}
\author{
Angelo Mazza, Antonio Punzo, Brian McGuire

Maintainer: Brian McGuire <mcguirebc@gmail.com>

}
\references{
 Mazza A, Punzo A, McGuire B. (2014). 
\pkg{KernSmoothIRT}: An {R} Package for Kernel Smoothing in Item Response Theory.
\emph{Journal of Statistical Software}, 58 6, 1-34. URL: \url{http://www.jstatsoft.org/v58/i06/}.

Ramsay, J.O. (2000). 
\pkg{TestGraf}: A program for the graphical analysis of multiple choice test and questionnaire data. 

Silverman, B.W. (1986). 
\emph{Density Estimation for Statistics and Data Analysis}

}