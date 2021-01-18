\name{residuals.dma}
\alias{residuals.dma}
\alias{residuals}

\title{Extracts Residuals from \code{dma} Model.}

\description{The function extracts residuals from the \code{\link{fDMA}} model.}

\usage{\method{residuals}{dma}(object, ...)}

\arguments{
\item{object}{an object of \code{dma} class}
\item{...}{not used}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dma")
r <- residuals(object=m1)
}
}
