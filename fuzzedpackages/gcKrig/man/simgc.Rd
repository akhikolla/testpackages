
\name{simgc}
\alias{simgc}
\title{Simulate Geostatistical Data from Gaussian Copula Model at Given Locations}
\usage{
simgc(locs, sim.n = 1, marginal, corr, longlat = FALSE)
}

\description{
	Simulate geostatistical data from Gaussian copula model at given locations.
	This function can simulate multiple datasets simultaneously.
}

\arguments{
\item{locs}{a numeric matrix or data frame of \emph{n-D} points with row denoting points.
First column is \emph{x} or longitude, second column is \emph{y} or latitude. The number of locations is
equal to the number of rows.}

\item{sim.n}{the number of simulation samples required.}

\item{marginal}{an object of class \code{\link{marginal.gc}} specifying the marginal distribution.}

\item{corr}{an object of class \code{\link{corr.gc}} specifying the correlation function.}

\item{longlat}{if FALSE, use Euclidean distance, if TRUE use great circle distance. Default is FALSE.}
}
\value{
A list of two elements:
        \item{data}{a numeric matrix with each row denoting a simulated data.}
        \item{locs}{the location of the simulated data, same as the input locs.}
}
\examples{
grid <- seq(0.05, 0.95, by = 0.1)
xloc <- expand.grid(x = grid, y = grid)[,1]
yloc <- expand.grid(x = grid, y = grid)[,2]
set.seed(12345)
sim1 <- simgc(locs = cbind(xloc,yloc), sim.n = 10, marginal = negbin.gc(mu = 5, od = 1),
              corr = matern.gc(range = 0.3, kappa = 0.5, nugget = 0.1))
#plot(sim1, index = 1)
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}

\keyword{Simulation}
