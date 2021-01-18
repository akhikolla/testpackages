\name{plotCounts}
\alias{plotCounts}
\title{Plot faecal egg count data}
\usage{
plotCounts(data, paired = TRUE, points = TRUE, 
    points.method = "jitter",  xlabel = "", 
    ylabel = "Faecal egg counts [epg]", ...)
}
\arguments{
\item{data}{a data.frame with two columns, the first column is before treatment counts, the second column is after treatment counts.}

\item{paired}{logical. If TRUE, uses the plot for the paired design. Otherwise uses the plot for the unpaired design.}

\item{points}{logical. If TRUE, add individual points for unpaired plot. Not used if \code{paired = TRUE}.}

\item{points.method}{string. It is used to separate coincident points if \code{points = TRUE}. One of \dQuote{jitter}, \dQuote{stack} or \dQuote{overplot}.}

\item{xlabel}{string. Label of x-axis.}

\item{ylabel}{String. Label of y-axis.}

\item{...}{Additional arguments for function \code{\link{xyplot}} if paired is TRUE, for function \code{\link{boxplot}} otherwise.}
}
\value{
A plot is drawn.
}
\description{
Plot egg count data to reflect changes between before and after treatment.
}
\author{
Craig Wang
}
\details{
For paired data, a xyplot is used. For unpiared data, a grouped boxplot is used.
}
\examples{
data(epgs)
plotCounts(epgs[,c("before","after")], paired = TRUE)
}
\keyword{hplot}
