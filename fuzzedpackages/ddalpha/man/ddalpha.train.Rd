\name{ddalpha.train}
\alias{ddalpha.train}

\alias{alpha}
\alias{polynomial}
\alias{knnlm}
\alias{maxD}

\alias{outsiders}

\title{
Train DD-Classifier
}
\description{
Trains the DD-classifier using a training sample according to given parameters. 
The DD-classifier is a non-parametric procedure that first transforms the training sample into the depth space calculating the depth of each point w.r.t each class (dimension of this space equals the number of classes in the training sample), and then constructs a separating rule in this depth space. 
If in the classification phase an object does not belong to the convex hull of at least one class (we mention such an object as an 'outsider'), it is mapped into the origin of the depth space and hence cannot be classified in the depth space. For these objects, after 'outsiderness' has been assured, an outsider treatment, i.e. a classification procedure functioning outside convex hulls of the classes is applied; it has to be trained too.

The current realization of the DD-classifier allows for several alternative outsider treatments; they involve different traditional classification methods, see 'Details' and 'Arguments' for parameters needed. 

The function allows for classification with \eqn{q\ge 2} classes, see \code{aggregation.method} in 'Arguments'.
}
\usage{
ddalpha.train(formula, data, subset,
              depth = "halfspace", 
              separator = "alpha", 
              outsider.methods = "LDA", 
              outsider.settings = NULL, 
              aggregation.method = "majority",
              pretransform = NULL,
              use.convex = FALSE,     
              seed = 0,
              ...)

}
\arguments{
  \item{formula}{
an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of the model. If not found in \code{data}, the variables of the model are taken from environment.
}
  \item{data}{
Matrix or data.frame containing training sample where each of \eqn{n} rows is one object of the training sample where first \eqn{d} entries are inputs and the last entry is output (class label).

A pre-calculated DD-plot may be used as \code{data} with \code{depth="ddplot"}.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in training the classifier.
}
  \item{depth}{
Character string determining which depth notion to use; the default value is \code{"halfspace"}. The list of the supported depths is given in section \emph{\bold{\code{Depths}}}. To use a custom depth, see topic \code{\link{Custom Methods}}. To use an outsider treatment only set \code{depth = NULL}.
}
  \item{separator}{
The method used for separation on the DD-plot; can be \code{"alpha"} (the default), \code{"polynomial"}, \code{"knnlm"} or \code{"maxD"}. See section \emph{\bold{\code{Separators}}} for the description of the separators and additional parameters. To use a custom separator, see topic \code{\link{Custom Methods}}.
}
  \item{outsider.methods}{
Vector of character strings each being a name of a basic outsider method for eventual classification; possible names are: \code{"LDA"} (the default), \code{"QDA"}, \code{"kNN"}, \code{"kNNAff"}, \code{"depth.Mahalanobis"}, \code{"RandProp"}, \code{"RandEqual"} and \code{"Ignore"}. Each method can be specified only once, replications are ignored. By specifying treatments in such a way only a basic treatment method can be chosen (by the name), and the default settings for each of the methods are applied, see 'Details'.
}
  \item{outsider.settings}{
List containing outsider treatments each described by a list of parameters including a name, see 'Details' and 'Examples'. Each method can be used multiply with (not necessarily) different parameters, just the name should be unique, entries with the repeating names are ignored.
}
  \item{aggregation.method}{
Character string determining which method to apply to aggregate binary classification results during multiclass classification; can be \code{"majority"} (the default) or \code{"sequent"}. If \code{"majority"}, \eqn{q(q-1)/2} (with \eqn{q} being the number of classes in the training sample) binary classifiers are trained, the classification results are aggregated using the majority voting, where classes with larger proportions in the training sample (eventually with the earlier entries in the \code{data}) are preferred when tied. If \code{"sequent"}, \eqn{q} binary 'one against all'-classifiers are trained and ties during the classification are resolved as before.
}
  \item{pretransform}{
indicates if the data has to be scaled before the learning procedure.
If the used depth method is affine-invariant and pretransform doesn't influence the result, the data won't be transformed (the parameter is ignored).

\describe{
  \item{NULL}{
applies no transformation to the data
  }
  \item{"1Mom", "1MCD"}{
the data is transformed with the common covariance matrix of the whole data
  }
  \item{"NMom", "NMCD"}{
the data is transformed w.r.t. each class using its covariance martix. The depths w.r.t. each class are calculated using the transformed data.
  }
  for the values \code{"1MCD", "NMCD"} \code{\link{covMcd}} is used to calculate the covariance matrix, and the parameter \code{mah.parMcd} is used.
  }
  
}
  \item{use.convex}{
Logical variable indicating whether outsiders should be determined exactly, i.e. as the points not contained in any of the convex hulls of the classes from the training sample (\code{TRUE}), or those having zero depth w.r.t. each class from the training sample (\code{FALSE}). For \code{depth =} \code{"zonoid"} both values give the same result.
}
  \item{seed}{
the random seed. The default value \code{seed=0} makes no changes.
}
  \item{...}{
The parameters for the depth calculating and separation methods.
}
}

\details{

\subsection{Depths}{

For \code{depth="ddplot"} the pre-calculated DD-plot shall be passed as \code{data}.

To use a custom depth, see topic \code{\link{Custom Methods}}.

To use an outsider treatment only set \code{depth = NULL}.

The following depths are supported:

\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

\code{\link{depth.zonoid}} for calculation of zonoid depth.

The additional parameters are described in the corresponding topics.

}

\subsection{Separators}{

The separators classify data on the 2-dimensional space of a DD-plot built using the depths.

To use a custom separator, see topic \code{\link{Custom Methods}}.

\subsection{alpha}{

Trains the DD\eqn{\alpha}-classifier (Lange, Mosler and Mozharovskyi, 2014; Mozharovskyi, Mosler and Lange, 2015). The DD\eqn{\alpha}-classifier constructs a linear separating rule in the polynomial extension of the depth space with the \eqn{\alpha}-procedure (Vasil'ev, 2003); maximum degree of the polynomial products is determined via cross-validation (in the depth space).

The additional parameters:
\describe{
  \item{max.degree}{
Maximum of the range of degrees of the polynomial depth space extension over which the \eqn{\alpha}-procedure is to be cross-validated; can be 1, 2 or 3 (default).
}
  \item{num.chunks}{
Number of chunks to split data into when cross-validating the \eqn{\alpha}-procedure; should be \eqn{>0}, and smaller than the total number of points in the two smallest classes when \code{aggregation.method =} \code{"majority"} and smaller than the total number of points in the training sample when \code{aggregation.method =} \code{"sequent"}. The default value is 10.
}
}
}

\subsection{polynomial}{

Trains the polynomial DD-classifier (Li, Cuesta-Albertos and Liu, 2012). The DD-classifier constructs a polynomial separating rule in the depth space; the degree of the polynomial is determined via cross-validation (in the depth space).

The additional parameters:
\describe{
  \item{max.degree}{
Maximum of the range of degrees of the polynomial over which the separator is to be cross-validated; can be in [1:10], the default value is 3.
}
  \item{num.chunks}{
Number of chunks to split data into when cross-validating the separator; should be \eqn{>0}, and smaller than the total number of points in the two smallest classes when \code{aggregation.method =} \code{"majority"} and smaller than the total number of points in the training sample when \code{aggregation.method =} \code{"sequent"}. The default value is 10.
}
}
}

\subsection{knnlm}{

Trains the \code{k}-nearest neighbours classifier in the depth space.

The additional parameters:
\describe{
  \item{knnrange}{
The maximal number of neighbours for kNN separation. The value is bounded by \eqn{2} and \eqn{n/2}.

\code{NULL} for the default value \eqn{10*(n^{1/q})+1}, where \eqn{n} is the number of objects, \eqn{q} is the number of classes. 

\code{"MAX"} for the maximum value \eqn{n/2}
}
}
}

\subsection{maxD}{
The \code{maximum depth} separator classifies an object to the class that provides it the largest depth value.
}

}

\subsection{Outsider treatment}{

An outsider treatment is a supplementary classifier for data that lie outside the convex hulls of all \eqn{q} training classes.
Available methods are: Linear Discriminant Analysis (referred to as "LDA"), see \code{\link{lda}}; \eqn{k}-Nearest-Neighbor Classifier ("kNN"), see \code{\link{knn}}, \code{\link{knn.cv}}; Affine-Invariant kNN ("kNNAff"), an affine-invariant version of the kNN, suited only for binary classification (some aggregation is used with multiple classes) and not accounting for ties (at all), but very fast by that; Maximum Mahalanobis Depth Classifier ("depth.Mahalanobis"), the outsider is referred to a class w.r.t. which it has the highest depth value scaled by (approximated) priors; Proportional Randomization ("RandProp"), the outsider is referred to a class randomly with probability equal to it (approximated) prior; Equal Randomization ("RandEqual"), the outsider is referred to a class randomly, chances for each class are equal; Ignoring ("Ignore"), the outsider is not classified, the string "Ignored" is returned instead.

An outsider treatment is specified by a list containing a name and parameters:

\code{name} is a character string, name of the outsider treatment to be freely specified; should be unique; is obligatory.

\code{method} is a character string, name of the method to use, can be \code{"LDA"}, \code{"kNN"}, \code{"kNNAff"}, \code{"depth.Mahalanobis"}, \code{"RandProp"}, \code{"RandEqual"} and \code{"Ignore"}; is obligatory.

\code{priors} is a numerical vector specifying prior probabilities of classes; class portions in the training sample are used by the default. \code{priors} is used in methods "LDA", "depth.Mahalanobis" and "RandProp".

\code{knn.k} is the number of the nearest neighbors taken into account; can be between \eqn{1} and the number of points in the training sample. Set to \eqn{-1} (the default) to be determined by the leave-one-out cross-validation. \code{knn.k} is used in method "kNN".

\code{knn.range} is the upper bound on the range over which the leave-one-out cross-validation is performed (the lower bound is \eqn{1}); can be between \eqn{2} and the number of points in the training sample \eqn{-1}. Set to \eqn{-1} (the default) to be calculated automatically accounting for number of points and dimension. \code{knn.range} is used in method "kNN".

\code{knnAff.methodAggregation} is a character string specifying the aggregation technique for method "kNNAff"; works in the same way as the function argument \code{aggregation.method}. \code{knnAff.methodAggregation} is used in method "kNNAff".

\code{knnAff.k} is the number of the nearest neighbors taken into account; should be at least \eqn{1} and up to the number of points in the training sample when \code{knnAff.methodAggregation =} \code{"sequent"}, and up to the total number of points in the training sample when \code{knnAff.methodAggregation =} \code{"majority"}. Set to \eqn{-1} (the default) to be determined by the leave-one-out cross-validation. \code{knnAff.k} is used in method "kNNAff".

\code{knnAff.range} is the upper bound on the range over which the leave-one-out cross-validation is performed (the lower bound is \eqn{1}); should be \eqn{>1} and smaller than the total number of points in the two smallest classes when \code{knnAff.methodAggregation =} \code{"majority"}, and \eqn{>1} and smaller than the total number of points in the training sample when \code{knnAff.methodAggregation =} \code{"sequent"}. Set to \eqn{-1} to be calculated automatically accounting for number of points and dimension. \code{knnAff.range} is used in method "kNNAff".

\code{mah.estimate} is a character string specifying which estimates to use when calculating the Mahalanobis depth; can be \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. \code{mah.estimate} is used in method "depth.Mahalanobis".

\code{mcd.alpha} is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used in method "depth.Mahalanobis" when \code{mah.estimate =} \code{"MCD"}.
}
}
\value{
Trained DD\eqn{\alpha}-classifier containing following - rather informative - fields:
\item{num.points}{Total number of points in the training sample.}
\item{dimension}{Dimension of the original space.}
\item{depth}{Character string determining which depth notion to use.}
\item{methodAggregation}{Character string determining which method to apply to aggregate binary classification results.}
\item{num.chunks}{Number of chunks data has been split into when cross-validating the \eqn{\alpha}-procedure.}
\item{num.directions}{Number of directions used for approximating the Tukey depth (when it is used).}
\item{use.convex}{Logical variable indicating whether outsiders should be determined exactly when classifying.}
\item{max.degree}{Maximum of the range of degrees of the polynomial depth space extension over which the \eqn{\alpha}-procedure has been cross-validated.}
\item{patterns}{Classes of the training sample.}
\item{num.classifiers}{Number of binary classifiers trained.}
\item{outsider.methods}{Treatments to be used to classify outsiders.}
}
\references{
Dyckerhoff, R., Koshevoy, G., and Mosler, K. (1996). Zonoid data depth: theory and computation. In: Prat A. (ed), \emph{COMPSTAT 1996. Proceedings in computational statistics}, Physica-Verlag (Heidelberg), 235--240.

Lange, T., Mosler, K., and Mozharovskyi, P. (2014). Fast nonparametric classification based on data depth. \emph{Statistical Papers} \bold{55} 49--69.

Li, J., Cuesta-Albertos, J.A., and Liu, R.Y. (2012). DD-classifier: Nonparametric classification procedure based on DD-plot. \emph{Journal of the American Statistical Association} \bold{107} 737--753.

Mozharovskyi, P. (2015). \emph{Contributions to Depth-based Classification and Computation of the Tukey Depth}. Verlag Dr. Kovac (Hamburg).

Mozharovskyi, P., Mosler, K., and Lange, T. (2015). Classifying real-world data with the DD\eqn{\alpha}-procedure. \emph{Advances in Data Analysis and Classification} \bold{9} 287--314.

Vasil'ev, V.I. (2003). The reduction principle in problems of revealing regularities I. \emph{Cybernetics and Systems Analysis} \bold{39} 686--694.
}
\seealso{
\code{\link{ddalpha.classify}} for classification using DD-classifier, 
\code{\link{depth.}} for calculation of depths, 
\code{\link{depth.space.}} for calculation of depth spaces, 
\code{\link{is.in.convex}} to check whether a point is not an outsider.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 200 training objects and 200 to test with
class1 <- mvrnorm(200, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(200, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
trainIndices <- c(1:100)
testIndices <- c(101:200)
propertyVars <- c(1:2)
classVar <- 3
trainData <- rbind(cbind(class1[trainIndices,], rep(1, 100)), 
                   cbind(class2[trainIndices,], rep(2, 100)))
testData <- rbind(cbind(class1[testIndices,], rep(1, 100)), 
                  cbind(class2[testIndices,], rep(2, 100)))
data <- list(train = trainData, test = testData)

# Train 1st DDalpha-classifier (default settings) 
# and get the classification error rate
ddalpha1 <- ddalpha.train(data$train)
classes1 <- ddalpha.classify(ddalpha1, data$test[,propertyVars])
cat("1. Classification error rate (defaults): ", 
    sum(unlist(classes1) != data$test[,classVar])/200, ".\n", sep = "")

# Train 2nd DDalpha-classifier (zonoid depth, maximum Mahalanobis 
# depth classifier with defaults as outsider treatment) 
# and get the classification error rate
ddalpha2 <- ddalpha.train(data$train, depth = "zonoid", 
                          outsider.methods = "depth.Mahalanobis")
classes2 <- ddalpha.classify(ddalpha2, data$test[,propertyVars], 
                               outsider.method = "depth.Mahalanobis")
cat("2. Classification error rate (depth.Mahalanobis): ", 
    sum(unlist(classes2) != data$test[,classVar])/200, ".\n", sep = "")

# Train 3rd DDalpha-classifier (100 random directions for the Tukey depth, 
# adjusted maximum Mahalanobis depth classifier 
# and equal randomization as outsider treatments) 
# and get the classification error rates
treatments <- list(list(name = "mahd1", method = "depth.Mahalanobis", 
                        mah.estimate = "MCD", mcd.alpha = 0.75, priors = c(1, 1)/2), 
                   list(name = "rand1", method = "RandEqual"))
ddalpha3 <- ddalpha.train(data$train, outsider.settings = treatments, 
                          num.direction = 100)
classes31 <- ddalpha.classify(ddalpha3, data$test[,propertyVars], 
                              outsider.method = "mahd1")
classes32 <- ddalpha.classify(ddalpha3, data$test[,propertyVars], 
                              outsider.method = "rand1")
cat("3. Classification error rate (by treatments):\n")
cat("   Error (mahd1): ", 
    sum(unlist(classes31) != data$test[,classVar])/200, ".\n", sep = "")
cat("   Error (rand1): ", 
    sum(unlist(classes32) != data$test[,classVar])/200, ".\n", sep = "")
    
# Train using some weird formula
ddalpha = ddalpha.train(
    I(mpg >= 19.2) ~ log(disp) + I(disp^2) + disp + I(disp * drat),
    data = mtcars, subset = (carb!=1), 
    depth = "Mahalanobis", separator = "alpha")
print(ddalpha) # make sure that the resulting table is what you wanted
CC = ddalpha.classify(ddalpha, mtcars)
sum((mtcars$mpg>=19.2)!= unlist(CC))/nrow(mtcars) # error rate
    
#Use the pre-calculated DD-plot
data = cbind(rbind(mvrnorm(n = 50, mu = c(0,0), Sigma = diag(2)),
                   mvrnorm(n = 50, mu = c(5,10), Sigma = diag(2)),
                   mvrnorm(n = 50, mu = c(10,0), Sigma = diag(2))),
             rep(c(1,2,3), each = 50))
plot(data[,1:2], col = (data[,3]+1))

ddplot = depth.space.Mahalanobis(data = data[,1:2], cardinalities = c(50,50,50))
ddplot = cbind(ddplot, data[,3])
ddalphaD = ddalpha.train(data = ddplot, depth = "ddplot", separator = "alpha")
c = ddalpha.classify(ddalphaD, ddplot[,1:3])
errors = sum(unlist(c) != data[,3])/nrow(data)
print(paste("Error rate: ",errors))

ddalpha = ddalpha.train(data = data, depth = "Mahalanobis", separator = "alpha")
c = ddalpha.classify(ddalpha, data[,1:2])
errors = sum(unlist(c) != data[,3])/nrow(data)
print(paste("Error rate: ",errors))
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }
