\name{gNormalize}
\alias{gNormalize}

\title{Normalizes a Numeric Matrix by Rows.}

\description{
For example, Google Trends data are given as numbers between 0 and 100. If the Users divide them by 100, they can be interpreted in a certain sense as probabilities. 

However, if there are such probabilities for several variables, sometimes it might be desirable to have the sum of these probabilities for all variables to sum up to 1. This function does not divide the values of an argument by 100, but rescales every row to sum up to 1. In other words, values in each row of an argument are divided by the sum of all values in this row. 
}
\usage{
gNormalize(data)
}

\arguments{
\item{data}{\code{\link[base]{matrix}}, observations are put in rows, and variables are grouped by columns}
}

\value{
\code{\link[base]{matrix}}
}

\examples{
gt <- gNormalize(trends)

gNormalize(rbind(c(0,1,2),c(1,2,3)))
}

\references{
Koop, G., Onorante, L., 2014. Macroeconomic nowcasting using Google probabilities. \url{https://goo.gl/ATsBN9}
}
