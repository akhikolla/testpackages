\name{print.grid.tvp}
\alias{print.grid.tvp}

\title{Prints \code{grid.tvp} Object.}

\description{The function prints selected outcomes obtained from \code{\link{grid.tvp}}.}

\usage{\method{print}{grid.tvp}(x, ...)}

\arguments{
\item{x}{an object of \code{grid.tvp} class}
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
grl <- c(0.99,0.98,0.97,0.96,0.95)
g <- grid.tvp(y=ld.wti,x=ld.drivers,V=1,grid.lambda=grl)

print(g)
}
}
