\name{summary.grid.dma}
\alias{summary.grid.dma}

\title{Summarizes Outcomes from \code{grid.dma} Objects.}

\description{The function summarizes outcomes obtained from \code{\link{grid.DMA}}.}

\usage{\method{summary}{grid.dma}(object, ...)}

\arguments{
\item{object}{an object of \code{grid.dma} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{\link{print.grid.dma}}.

Additionally, it finds the indices for a model minimizing Root Mean Squared Error (RMSE) and for a model minimizing Mean Absolute Error (MAE). 
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
gra <- c(0.99,0.98,0.97)
grl <- c(0.99,0.95)
g1 <- grid.DMA(y=ld.wti,x=ld.drivers,grid.alpha=gra,grid.lambda=grl,initvar=1)
g2 <- grid.DMA(y=ld.wti,x=ld.drivers,grid.alpha=gra,grid.lambda=grl,initvar=1,model="dms")
}
\dontshow{
gra <- c(0.99,0.98,0.97)
grl <- c(0.99,0.95)
s <- ld.drivers[1:20,1:4]
g1 <- grid.DMA(y=ld.wti[1:20,1],x=s,grid.alpha=gra,grid.lambda=grl,initvar=1)
g2 <- grid.DMA(y=ld.wti[1:20,1],x=s,grid.alpha=gra,grid.lambda=grl,initvar=1,model="dms")
}
summary(g1)
summary(g2)
}
