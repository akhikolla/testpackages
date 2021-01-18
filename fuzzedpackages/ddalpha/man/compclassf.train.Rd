\name{compclassf.train}
\alias{compclassf.train}
\title{
Functional Componentwise Classifier
}
\description{
Trains the functional componentwise classifier
}
\usage{
compclassf.train (dataf, labels, subset,
                  to.equalize = TRUE, 
                  to.reduce = TRUE, 
                  classifier.type = c("ddalpha", "maxdepth", "knnaff", "lda", "qda"), 
                  ...)
}

\arguments{
  \item{dataf}{
list containing lists (functions) of two vectors of equal length, named "args" and "vals": arguments sorted in ascending order and corresponding them values respectively
}
  \item{labels}{
list of output labels of the functional observations
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in training the classifier.
}
  \item{to.equalize}{
Adjust the data to have equal (the largest) argument interval.
}
  \item{to.reduce}{
If the data spans a subspace only, project on it (by PCA).
}
  \item{classifier.type}{
the classifier which is used on the transformed space. The default value is 'ddalpha'.
}
  \item{\dots}{
additional parameters, passed to the classifier, selected with parameter \code{classifier.type}.
}
}
\details{
The finite-dimensional space is directly constructed from the observed values. 
Delaigle, Hall and Bathia (2012) consider (almost) all sets of discretization points
that have a given cardinality. 

The usual classifiers are then trained on the constructed finite-dimensional space.
}
\value{
Trained functional componentwise classifier
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Delaigle, A., Hall, P., and Bathia, N. (2012). Componentwise classification and clustering of functional data. \emph{Biometrika} \bold{99} 299--313.
}

\seealso{
 \code{\link{compclassf.classify}} for classification using functional componentwise classifier,
 
 \code{\link{ddalphaf.train}} to train the functional DD-classifier,
 
 \code{\link{dataf.*}} for functional data sets included in the package.
 
}
\examples{

\dontrun{
## load the Growth dataset
dataf = dataf.growth()

learn = c(head(dataf$dataf, 49), tail(dataf$dataf, 34))
labels =c(head(dataf$labels, 49), tail(dataf$labels, 34)) 
test = tail(head(dataf$dataf, 59), 10)    # elements 50:59. 5 girls, 5 boys

c = compclassf.train (learn, labels, classifier.type = "ddalpha")

classified = compclassf.classify(c, test)

print(unlist(classified))

}
}

\keyword{ functional }
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }
