\name{grid.DMA}
\alias{grid.DMA}

\title{Computes \code{\link{fDMA}} Function for Multiple Values of \code{alpha} and \code{lambda}.}

\description{
Sometimes it is necessary to consider various values of parameters \code{alpha} and \code{lambda} in Dynamic Model Averaging (or Dynamic Model Selection, etc.). This function computes \code{\link{fDMA}} function for all combinations of \code{alpha} and \code{lambda} for given grids.

This function is a wrapper of \code{\link{fDMA}}.
}

\usage{
grid.DMA(y,x,grid.alpha,grid.lambda,initvar,W=NULL,initial.period=NULL,V.meth=NULL,
kappa=NULL,gprob=NULL,omega=NULL,model=NULL,parallel.grid=NULL,m.prior=NULL,
mods.incl=NULL,DOW=NULL,DOW.nmods=NULL,DOW.type=NULL,DOW.limit.nmods=NULL,
forced.models=NULL,forbidden.models=NULL,forced.variables=NULL,bm=NULL,
small.c=NULL,av=NULL)
}

\arguments{
\item{y}{see \code{\link{fDMA}}}
\item{x}{see \code{\link{fDMA}}}
\item{grid.alpha}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of different values of \code{alpha}}
\item{grid.lambda}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of different values of \code{lambda} or a \code{\link[base]{list}} of numeric vectors for multiple \code{lambda} in one model (see \code{\link{fDMA}})}
\item{initvar}{see \code{\link{fDMA}}}
\item{W}{see \code{\link{fDMA}}}
\item{initial.period}{see \code{\link{fDMA}}}
\item{V.meth}{see \code{\link{fDMA}}}
\item{kappa}{see \code{\link{fDMA}}}
\item{gprob}{see \code{\link{fDMA}}}
\item{omega}{see \code{\link{fDMA}}}
\item{model}{see \code{\link{fDMA}}}
\item{parallel.grid}{optional, \code{\link[base]{logical}}, indicate whether parallel computations should be used, by default \code{parallel.grid=FALSE}}
\item{m.prior}{see \code{\link{fDMA}}}
\item{mods.incl}{see \code{\link{fDMA}}}
\item{DOW}{see \code{\link{fDMA}}}
\item{DOW.nmods}{see \code{\link{fDMA}}}
\item{DOW.type}{see \code{\link{fDMA}}}
\item{DOW.limit.nmods}{see \code{\link{fDMA}}}
\item{forced.models}{see \code{\link{fDMA}}}
\item{forbidden.models}{see \code{\link{fDMA}}}
\item{forced.variables}{see \code{\link{fDMA}}}
\item{bm}{see \code{\link{fDMA}}}
\item{small.c}{see \code{\link{fDMA}}}
\item{av}{see \code{\link{fDMA}}}
}

\value{an object of class \code{grid.dma}, \code{\link[base]{list}} of

\item{$models}{\code{\link[base]{list}} of \code{\link[base]{list}} of models}
\item{$RMSE}{\code{\link[base]{matrix}} with Root Mean Squared Error (RMSE) for all estimated models}
\item{$MAE}{\code{\link[base]{matrix}} with Mean Absolute Error (MAE) for all estimated models}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
gra <- c(0.99,0.98,0.97)
grl <- c(0.99,0.95)
g1 <- grid.DMA(y=ld.wti,x=ld.drivers,grid.alpha=gra,grid.lambda=grl,initvar=10)
}
\dontshow{
gra <- c(0.99,0.98,0.97)
grl <- c(0.99,0.95)
g <- grid.DMA(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],grid.alpha=gra,grid.lambda=grl,initvar=10)
}
# extract model with alpha=0.97 and lambda=0.95
model <- g$models[[3]][[2]]

# models with various multiple lambdas
\dontrun{
gra <- c(0.99,0.98,0.97)
grl <- list(c(0.99,0.95,0.90),c(0.99,0.98,0.97,0.96,0.95))
g2 <- grid.DMA(y=ld.wti,x=ld.drivers,grid.alpha=gra,grid.lambda=grl,initvar=10)
}
}

\seealso{\code{\link{print.grid.dma}}, \code{\link{summary.grid.dma}}, \code{\link{plot.grid.dma}}.}
