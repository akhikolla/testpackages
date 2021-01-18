\name{spnn-package}
\alias{spnn-package}
\alias{spnn}
\docType{package}
\title{
\packageTitle{spnn}
}
\description{
\packageDescription{spnn}
}
\details{
The package exports 4 main functions:

\itemize{
\item{\link{spnn.learn}}{ Create or update a Scale Invariant Probabilistic Neural Network. }
\item{\link{spnn.predict}}{ Estimates the category probabilities of new observations using a fitted SPNN. }
\item{\link{cspnn.learn}}{ Create or update a Condensed Scale Invariant Probabilistic Neural Network. }
\item{\link{cspnn.predict}}{ Estimates the category probabilities of new observations using a fitted CSPNN. }
}
}
\author{
\packageMaintainer{spnn}
}
\references{
[1] Specht, Donald F. "Probabilistic neural networks." Neural networks 3.1 (1990): 109-118.

[2] Specht, Donald F. "Enhancements to probabilistic neural networks." Neural Networks, 1992.IJCNN., International Joint Conference on. Vol. 1. IEEE, 1992.

[3] Ebrahimi, Romin "Scale Invariant Probabilistic Neural Networks." The University of Texas, 2018 \url{https://repositories.lib.utexas.edu/handle/2152/65166}
}
\keyword{spnn}
\keyword{neural}
\keyword{network}
\keyword{probability}
\keyword{pnn}
\seealso{
\code{\link{spnn.learn}},
\code{\link{spnn.predict}},
\code{\link{cspnn.learn}},
\code{\link{cspnn.predict}}
}
\examples{
library(spnn)
library(datasets)

data(iris)

# shuffle the iris data set
indexRandom <- sample(1:nrow(iris), size = nrow(iris), replace = FALSE)

# use 100 observations for training set
trainData <- iris[indexRandom[1:100],]

# use remaining observations for testing
testData <- iris[indexRandom[101:length(indexRandom)],]

# fit spnn
spnn <- spnn.learn(set = trainData, category.column = 5)

# estimate probabilities
predictions <- spnn.predict(nn = spnn, newData = testData[,1:4])

# reference matrix must be supplied
# this is not the optimal reference matrix
# this matrix is provided as a simple example
xr <- matrix(c(c(5.00, 3.41, 1.44, 0.24),
               c(5.88, 2.75, 4.23, 1.30),
               c(6.61, 2.97, 5.59, 2.01)),
             nrow = length(unique(trainData$Species)),
             ncol = ncol(trainData) - 1,
             byrow = TRUE)

# fit cspnn
cspnn <- cspnn.learn(set = trainData, xr = xr, category.column = 5)

# estimate probabilities
predictions <- cspnn.predict(nn = cspnn, newData = testData[,1:4])

}
