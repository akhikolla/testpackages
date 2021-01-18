\name{onevar}
\alias{onevar}

\title{Creates a \code{\link[base]{matrix}} of one-variable models.}

\description{
This function simplifies working with one-variable models in, for example, \code{\link{fDMA}}. It produces a \code{\link[base]{matrix}} corresponding to the set of models consisting of models with a constant and just one extra variable, and a model with a constant only. 
}

\usage{
onevar(x)
}

\arguments{
\item{x}{\code{\link[base]{matrix}} of independent variables, see \code{mods.incl} in \code{\link{fDMA}}}
}

\value{\code{\link[base]{matrix}}, 
inclusion of a variable is indicated by 1, omitting by 0}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]

mds <- diag(1,ncol(ld.drivers),ncol(ld.drivers))
mds <- cbind(rep(1,ncol(ld.drivers)),mds)
mds <- rbind(rep(0,ncol(mds)),mds)
mds[1,1] <- 1

m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,mods.incl=mds)

# Equivalently: 

m2 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10,mods.incl=onevar(ld.drivers))
}
}
