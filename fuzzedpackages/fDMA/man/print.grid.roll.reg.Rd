\name{print.grid.roll.reg}
\alias{print.grid.roll.reg}

\title{Prints \code{grid.roll.reg} Object.}

\description{The function prints selected outcomes obtained from \code{\link{grid.roll.reg}}.}

\usage{\method{print}{grid.roll.reg}(x, ...)}

\arguments{
\item{x}{an object of \code{grid.roll.reg} class}
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
grw <- c(50,100,150)
g <- grid.roll.reg(y=ld.wti,x=ld.drivers,grid.window=grw)

print(g)
}
}
