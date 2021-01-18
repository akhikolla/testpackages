\name{reduce.size}
\alias{reduce.size}

\title{Reduces the Size of \code{\link{fDMA}} or \code{\link{grid.DMA}} Outcomes.}

\description{
This functions reduces the size of \code{dma} or \code{grid.dma} object. 
}

\details{The information corresponding to each sub-model is erased. In particular, for the object produced by \code{\link{fDMA}} \code{$models} is reduced to one-row matrix to keep only  \code{\link[base]{colnames}}, and \code{$postmod}, \code{$yhat.all.mods} and \code{$p.dens.} are replaced by \code{\link[base]{NA}}. It can be useful if large number of models is considered.}

\usage{
reduce.size(dma.object)
}

\arguments{
\item{dma.object}{\code{dma} or \code{grid.dma} object}
}

\value{\code{dma} or \code{grid.dma} object, 
with the information corresponding to each sub-model erased}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10)
m2 <- reduce.size(m1)}
}

\seealso{\code{\link{fDMA}}, \code{\link{grid.DMA}}.}
