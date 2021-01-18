\name{summary.grid.tvp}
\alias{summary.grid.tvp}

\title{Summarizes Outcomes from \code{grid.tvp} Objects.}

\description{The function summarizes outcomes obtained from \code{\link{grid.tvp}}.}

\usage{\method{summary}{grid.tvp}(object, ...)}

\arguments{
\item{object}{an object of \code{grid.tvp} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{\link{print.grid.tvp}}.

Additionally, it finds the model minimizing Root Mean Squared Error (RMSE) and minimizing Mean Absolute Error (MAE). 
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
grl <- c(0.99,0.98,0.97,0.96,0.95)
g <- grid.tvp(y=ld.wti,x=ld.drivers,V=1,grid.lambda=grl)

summary(g)
}
}

