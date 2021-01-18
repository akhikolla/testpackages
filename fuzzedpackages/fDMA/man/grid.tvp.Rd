\name{grid.tvp}
\alias{grid.tvp}

\title{Computes \code{\link{tvp}} Function for Multiple Values of \code{lambda}.}

\description{
Sometimes it is necessary to consider various values of parameter \code{lambda} in Time-Varying Parameters Regression. This function computes \code{\link{tvp}} function for all values of \code{lambda} for a given grid.

This function is a wrapper of \code{\link{tvp}}.
}

\usage{
grid.tvp(y,x,V,grid.lambda,W=NULL,kappa=NULL,parallel.grid=NULL,c=NULL)
}

\arguments{
\item{y}{see \code{\link{tvp}}}
\item{x}{see \code{\link{tvp}}}
\item{V}{see \code{\link{tvp}}}
\item{grid.lambda}{a numeric \code{\link[base]{vector}} of different values of \code{lambda}, see \code{\link{tvp}}}
\item{W}{optional, see \code{\link{tvp}}}
\item{kappa}{optional, see \code{\link{tvp}}}
\item{parallel.grid}{optional, \code{\link[base]{logical}}, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
\item{c}{optional, see \code{\link{tvp}}}
}

\value{an object of class \code{grid.tvp}, \code{\link[base]{list}} of

\item{$models}{\code{\link[base]{list}} of \code{\link{tvp}} objects}
\item{$fq}{\code{\link[base]{matrix}} with Root Mean Squared Error (RMSE)  and Mean Absolute Error (MAE) for all estimated models}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
grl <- c(0.99,0.98,0.97,0.96,0.95)
g <- grid.tvp(y=ld.wti,x=ld.drivers,V=1,grid.lambda=grl)

# extract model with lambda=0.95
model <- g$models[[5]]
}
}

\seealso{\code{\link{print.grid.tvp}}, \code{\link{summary.grid.tvp}}, \code{\link{plot.grid.tvp}}.}
