\name{displace}
\alias{displace}
\title{Estimate horizontal displacement and velocity between two dates}
\description{
Returns an \R \code{data.frame} containing the crustal displacement in meters and velocity in mm/yr at a specified location from time t1 to time t2. The estimated displacement equals the velocity at this location multiplied by the time difference (t2 - t1) plus all coseismic and postseismic motion that has occurred between these two times.
}
\usage{
displace(xy, t0, t1, iopt)
}
\arguments{
  \item{xy}{A matrix of longitude and latitude coordinates in decimal degrees}
  \item{t0}{The starting date}
  \item{t1}{The ending date}
  \item{iopt}{The key identifying the reference frame for the input values. See \code{\link{iopt}} for a list of supported values.}
}
\value{
  A \code{data.frame} with 6 columns:
    \item{de}{Eastward displacement (meters)}
    \item{dn}{Northward displacement (meters)}
    \item{du}{Upward displacement (meters)}
    \item{ve}{Eastward velocity (mm/yr)}
    \item{vn}{Northward velocity (mm/yr)}
    \item{vu}{Upward velocity (mm/yr)}
}
\examples{
beta <- matrix(data=c(-121.7722, 36.6698), nrow=1)
t0 <- as.Date("1989-10-16")
t1 <- as.Date("1989-10-18")
displace(beta, t0, t1, 1)
}
\seealso{
  \code{\link{iopt}}
}
