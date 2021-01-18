\name{coef.dma}
\alias{coef.dma}
\alias{coef}

\title{Extracts Averaged Coefficients from \code{dma} Model.}

\description{The function extracts the expected values of regression coefficients from the \code{\link{fDMA}} model.}

\usage{\method{coef}{dma}(object, ...)}

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
c <- coef(object=m1)
}
}
