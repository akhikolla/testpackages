\name{spnn.predict}
\alias{spnn.predict}
\title{
spnn.predict
}
\description{

Estimates the category probabilities of new observations using a fitted SPNN.
}
\usage{
spnn.predict(nn, newData)
}

\arguments{
  \item{nn}{

A trained Scaled Invariant Probabilistic Neural Network.
}
  \item{newData}{
A matrix of new observations where each row represents a single observation vector.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Given a trained Scale Invariant Probabilistic Neural Network and new data, the function spnn.predict returns the category with the highest probability and the probability estimates for each category.
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
\code{\link{spnn-package}}, \code{\link{spnn.learn}}, \code{\link{iris}}
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

}
