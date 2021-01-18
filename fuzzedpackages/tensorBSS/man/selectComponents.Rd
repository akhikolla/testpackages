\name{selectComponents}
\alias{selectComponents}
\title{Select the Most Informative Components}
\description{
Takes an array of observations as an input and outputs a subset of the components having the most extreme kurtoses.
}
\usage{
selectComponents(x, first = 2, last = 2)
}
\arguments{
  \item{x}{Numeric array of an order at least two. It is assumed that the last dimension corresponds to the sampling units.}
  \item{first}{Number of components with maximal kurtosis to be selected. Can equal zero but the total number of components selected must be at least two.}
  \item{last}{Number of components with minimal kurtosis to be selected. Can equal zero but the total number of components selected must be at least two.}
}
\details{
In independent component analysis (ICA) the components having the most extreme kurtoses are often thought to be also the most informative. With this viewpoint in mind the function \code{selectComponents} selects from \code{x} \code{first} components having the highest kurtosis and \code{last} components having the lowest kurtoses and outputs them as a standard data matrix for further analysis.
}
\value{
Data matrix with rows corresponding to the observations and the columns correponding to the \code{first} + \code{last} selected components in decreasing order with respect to kurtosis. The names of the components in the output matrix correspond to the indices of the components in the original array \code{x}.
}
\author{
Joni Virta
}

\examples{
data(zip.train)
x <- zip.train

rows <- which(x[, 1] == 0 | x[, 1] == 1)
x0 <- x[rows, 2:257]

x0 <- t(x0)
dim(x0) <- c(16, 16, 2199)

tfobi <- tFOBI(x0)
comp <- selectComponents(tfobi$S)
head(comp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

