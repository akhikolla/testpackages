\name{predict.dma}
\alias{predict.dma}
\alias{predict}

\title{Computes Predictions from \code{dma} Model.}

\description{The function computes predictions based on the model obtained from \code{\link{fDMA}}.}

\usage{\method{predict}{dma}(object, newdata, type, ...)}

\arguments{
\item{object}{an object of \code{dma} class}
\item{newdata}{a \code{matrix} as \code{x} object in \code{\link{fDMA}}}
\item{type}{\code{type="backward"} computes predictions of \code{y} with the already estimated coefficients, but with \code{x} given by \code{newdata}, \code{type="forward"} computes predictions of \code{y} with the coefficients estimated in the last period, for various combinations of \code{x} given in rows of \code{newdata}}
\item{...}{not used}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=1,model="dma")
p1 <- predict(object=m1,newdata=ld.drivers,type="backward")
p2 <- predict(object=m1,newdata=ld.drivers[1,],type="forward")
p3 <- predict(object=m1,newdata=ld.drivers[1:3,],type="forward")
}
}
