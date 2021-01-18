\name{print.grid.dma}
\alias{print.grid.dma}

\title{Prints \code{grid.dma} Object.}

\description{The function prints selected outcomes obtained from \code{\link{grid.DMA}}.}

\usage{\method{print}{grid.dma}(x, ...)}

\arguments{
\item{x}{an object of \code{grid.dma} class}
\item{...}{not used}
}

\details{
The function prints Root Mean Squared Error (RMSE) and Mean Absolute Error (MAE) for all estimated models.
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
print(g1)
print(g2)
}
