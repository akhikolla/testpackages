\name{print.tvp}
\alias{print.tvp}

\title{Prints \code{tvp} Object.}

\description{The function prints selected outcomes obtained from \code{\link{tvp}}.}

\usage{\method{print}{tvp}(x, ...)}

\arguments{
\item{x}{an object of \code{tvp} class}
\item{...}{not used}
}

\details{
The function prints mean regression coefficients from the analyzed period, Root Mean Squared Error (RMSE) and Mean Absolute Error (MAE) from the estimated model. 
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
tvp <- tvp(y=ld.wti,x=ld.drivers,V=1,lambda=0.99)
print(tvp)
}
