\name{cspnn.predict}
\alias{cspnn.predict}
\title{
cspnn.predict
}
\description{

Estimates the category probabilities of new observations using a fitted CSPNN.
}
\usage{
cspnn.predict(nn, newData)
}

\arguments{
  \item{nn}{

A trained Condensed Scaled Invariant Probabilistic Neural Network.
}
  \item{newData}{
A matrix of new observations where each row represents a single observation vector.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Given a trained Condensed Scale Invariant Probabilistic Neural Network and new data, the function cspnn.predict returns the category with the highest probability and the probability estimates for each category.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of the guessed categories and the probability estimates of each category.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{spnn-package}}, \code{\link{cspnn.learn}}, \code{\link{iris}}
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
