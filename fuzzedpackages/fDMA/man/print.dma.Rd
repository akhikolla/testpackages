\name{print.dma}
\alias{print.dma}
\alias{print}

\title{Prints \code{dma} Object.}

\description{The function prints selected outcomes obtained from \code{\link{fDMA}}.}

\usage{\method{print}{dma}(x, ...)}

\arguments{
\item{x}{an object of \code{dma} class}
\item{...}{not used}
}

\details{
The function prints parameters of an argument \code{x}, Root Mean Squared Error (RMSE) and Mean Absolute Error (MAE) from the estimated model. It also shows the number of observations, the number of models in averaging (selecting) procedure and the number of variables (including constant) used in the model. The number of models does not include the increase, if multiple \code{lambda} is used. The function also shows forecast quality measures for alternative forecasting methods, i.e., naive forecast (see also \code{\link{altf}}) and, if computed, for Auto ARIMA \code{\link[forecast]{auto.arima}}.
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dma")
m2 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dms")
}
\dontshow{
m1 <- fDMA(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],alpha=0.99,lambda=0.99,initvar=1,model="dma")
m2 <- fDMA(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],alpha=0.99,lambda=0.99,initvar=1,model="dms")
}
print(m1)
print(m2)
}
