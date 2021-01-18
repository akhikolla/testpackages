\name{grid.roll.reg}
\alias{grid.roll.reg}

\title{Computes \code{\link{roll.reg}} Function for Multiple Values of \code{window}.}

\description{
Sometimes it is necessary to consider various values of parameter \code{window} in Rolling Regression. This function computes \code{\link{roll.reg}} function for all values of \code{window} for a given grid.

This function is a wrapper of \code{\link{roll.reg}}.
}

\usage{
grid.roll.reg(y,x=NULL,grid.window,parallel.grid=NULL,c=NULL)
}

\arguments{
\item{y}{see \code{\link{roll.reg}}}
\item{x}{see \code{\link{roll.reg}}}
\item{grid.window}{a numeric \code{\link[base]{vector}} of different values of \code{window}, see \code{\link{roll.reg}}}
\item{parallel.grid}{optional, \code{\link[base]{logical}}, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
\item{c}{optional, see \code{\link{roll.reg}}}
}

\value{an object of class \code{grid.roll.reg}, \code{\link[base]{list}} of

\item{$models}{\code{\link[base]{list}} of \code{reg} objects}
\item{$fq}{\code{\link[base]{matrix}} with Root Mean Squared Error (RMSE)  and Mean Absolute Error (MAE) for all estimated models}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
grw <- c(50,100,150)
g <- grid.roll.reg(y=ld.wti,x=ld.drivers,grid.window=grw)

# extract model with window=100
model <- g$models[[2]]
}
}

\seealso{\code{\link{print.grid.roll.reg}}, \code{\link{summary.grid.roll.reg}}, \code{\link{plot.grid.roll.reg}}.}
