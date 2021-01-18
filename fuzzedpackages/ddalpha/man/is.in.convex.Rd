\name{is.in.convex}
\alias{is.in.convex}
\title{
Check Outsiderness
}
\description{
Checks the belonging to at least one of class convex hulls of the training sample.
}
\usage{
is.in.convex(x, data, cardinalities, seed = 0)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose belonging to convex hulls is to be checked; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects, representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{seed}{
the random seed. The default value \code{seed=0} makes no changes.
}
}
\details{
Checks are conducted w.r.t. each separate class in \code{data} using the simplex algorithm, taken from the C++ implementation of the zonoid depth calculation  by Rainer Dyckerhoff.
}
\value{
Matrix of \code{number of objects} rows and \code{number of classes} columns, containing \code{1} if an object belongs to the convex hull of the corresponding class, and \code{0} otherwise.
}
\references{
Dyckerhoff, R., Koshevoy, G., and Mosler, K. (1996). Zonoid data depth: theory and computation. In: Prat A. (ed), \emph{COMPSTAT 1996. Proceedings in computational statistics}, Physica-Verlag (Heidelberg), 235--240.
}
\author{
Implementation of the simplex algorithm is taken from the algorithm for computation of zonoid depth (Dyckerhoff, Koshevoy and Mosler, 1996) that has been implemented in C++ by Rainer Dyckerhoff.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 400 training objects and 1000 to test with
class1 <- mvrnorm(700, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(700, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
trainIndices <- c(1:200)
testIndices <- c(201:700)
propertyVars <- c(1:2)
classVar <- 3
trainData <- rbind(cbind(class1[trainIndices,], rep(1, 200)), 
                   cbind(class2[trainIndices,], rep(2, 200)))
testData <- rbind(cbind(class1[testIndices,], rep(1, 500)), 
                  cbind(class2[testIndices,], rep(2, 500)))
data <- list(train = trainData, test = testData)

# Count outsiders
numOutsiders = sum(rowSums(is.in.convex(data$test[,propertyVars], 
                                data$train[,propertyVars], c(200, 200))) == 0)
cat(numOutsiders, "outsiders found in the testing sample.\n")
}
\keyword{ multivariate }
