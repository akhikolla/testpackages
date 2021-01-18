\name{rvi}
\alias{rvi}

\title{Extracts Relative Variable Importances from \code{\link{fDMA}} Model.}

\description{
This functions extracts posterior inclusion probabilities for independent variables from \code{dma} object. 
}

\usage{
rvi(dma.object)
}

\arguments{
\item{dma.object}{\code{dma} object}
}

\value{\code{matrix} of posterior inclusion probabilities for independent variables}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10)
r <- rvi(m1)}
}
