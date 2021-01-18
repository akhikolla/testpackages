\name{PartitionParameters}
\alias{PartitionParameters}
\title{Partition parameters}
\usage{
PartitionParameters(windowsweep, type = c("smooth", "flat")[1],
  clusterwidth = 1)
}
\arguments{
  \item{windowsweep}{a \code{windowsweep} object, i.e. the
  output of the \code{\link{WindowSweep}} function.}

  \item{type}{type of estimate to output, whether "smooth"
  or "flat".}

  \item{clusterwidth}{the temporal range within which
  changepoints are considered to be within the same cluster
  (for a "smooth" BCPA).}
}
\value{
a data frame containing the three estimates: \code{mu.hat},
\code{s.hat}, \code{rho.hat}.
}
\description{
Partitions - and, ultimately, estimates - all parameters of
a BCPA, either as a rolling average (smooth BCPA), or as
constant values within fixed change points (flat BCPA).
}
\author{
Eliezer Gurarie
}
\seealso{
used in \code{\link{ChangePointSummary}} and
\code{\link{PhasePlot}}
}

