\name{print.reg}
\alias{print.reg}

\title{Prints \code{reg} Object.}

\description{The function prints selected outcomes obtained from \code{\link{roll.reg}} and \code{\link{rec.reg}}.}

\usage{\method{print}{reg}(x, ...)}

\arguments{
\item{x}{an object of \code{reg} class}
\item{...}{not used}
}

\details{
The function prints mean regression coefficients from the analyzed period, Root Mean Squared Error (RMSE) and Mean Absolute Error (MAE) from the estimated model. For \code{\link{roll.reg}} it also shows the size of a rolling window. 
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
roll <- roll.reg(y=ld.wti,x=ld.drivers,window=100)
rec <-  rec.reg(y=ld.wti,x=ld.drivers)
print(roll)
print(rec)
}
