\name{ddalphaf.train}
\alias{ddalphaf.train}
\title{
Functional DD-Classifier
}
\description{
Trains the functional DD-classifier
}
\usage{
ddalphaf.train(dataf, labels, subset, 
                adc.args = list(instance = "avr", 
                               numFcn = -1, 
                               numDer = -1), 
                classifier.type = c("ddalpha", "maxdepth", "knnaff", "lda", "qda"), 
                cv.complete = FALSE, 
                maxNumIntervals = min(25, ceiling(length(dataf[[1]]$args)/2)),
                seed = 0,
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
  \item{adc.args}{
Represents a function sample as a multidimensional (dimension=\code{"numFcn"+"numDer"}) 
one averaging (\code{instance = "avr"}) or evaluating (\code{instance = "val"}) for that each function and it derivative on \code{"numFcn"} 
(resp. \code{"numDer"}) equal nonoverlapping covering intervals

First two named \code{"args"} and \code{"vals"} are arguments sorted in 
ascending order and having same bounds for all functions and 
corresponding them values respectively

  \describe{
    \item{instance}{
type of discretizing the functions: \cr
"avr" - by averaging over intervals of the same length \cr
"val" - by taking values on equally-spaced grid
    }
    \item{numFcn}{
  number of function intervals
    }  
    \item{numDer}{
  number of first-derivative intervals 
    }
  }
  Set \code{numFcn} and \code{numDer} to -1 to apply cross-validation.
  
  Set  \code{adc.args} to a list of "adc.args" objects to cross-validate only over these values.
}
  \item{classifier.type}{
the classifier which is used on the transformed space. The default value is 'ddalpha'.
}
  \item{cv.complete}{
T: apply complete cross-validation\cr
F: restrict cross-validation by Vapnik-Chervonenkis bound
}
  \item{maxNumIntervals}{
maximal number of intervals for cross-validation ( max(numFcn + numDer) = maxNumIntervals )
}
  \item{seed}{
the random seed. The default value \code{seed=0} makes no changes.
}
  \item{\dots}{
additional parameters, passed to the classifier, selected with parameter \code{classifier.type}.
}
}
\details{
The functional DD-classifier is fast nonparametric procedure for classifying functional data. It consists of
a two-step transformation of the original data plus a classifier operating on a low-dimensional
hypercube. The functional data are first mapped into a finite-dimensional location-slope space
and then transformed by a multivariate depth function into the DD-plot, which is a subset of
the unit hypercube. This transformation yields a new notion of depth for functional data. Three
alternative depth functions are employed for this, as well as two rules for the final classification. 
The resulting classifier is cross-validated over a small range of parameters only,
which is restricted by a Vapnik-Cervonenkis bound. The entire methodology does not involve
smoothing techniques, is completely nonparametric and allows to achieve Bayes optimality under
standard distributional settings. It is robust and efficiently computable.
}
\value{
Trained functional DD-classifier
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Mosler, K. and Mozharovskyi, P. (2017). Fast DD-classification of functional data. \emph{Statistical Papers} \bold{58} 1055--1089.

Mozharovskyi, P. (2015). \emph{Contributions to Depth-based Classification and Computation of the Tukey Depth}. Verlag Dr. Kovac (Hamburg).
}

\seealso{
 \code{\link{ddalphaf.classify}} for classification using functional DD\eqn{\alpha}-classifier,
 \code{\link{compclassf.train}} to train the functional componentwise classifier,
 \code{\link{dataf.*}} for functional data sets included in the package.
}
\examples{

\dontrun{

## load the Growth dataset
dataf = dataf.growth()

learn = c(head(dataf$dataf, 49), tail(dataf$dataf, 34))
labels= c(head(dataf$labels, 49), tail(dataf$labels, 34)) 
test  = tail(head(dataf$dataf, 59), 10)    # elements 50:59. 5 girls, 5 boys

#cross-validate over the whole variants up to dimension 3
c1 = ddalphaf.train (learn, labels, classifier.type = "ddalpha", maxNumIntervals = 3)

classified1 = ddalphaf.classify(c1, test)

print(unlist(classified1))
print(c1$adc.args) 

# cross-validate over these two variants
c2 = ddalphaf.train (learn, labels, classifier.type = "ddalpha", 
                     adc.args = list(
                       list(instance = "avr", 
                            numFcn = 1, 
                            numDer = 2),
                       list(instance = "avr", 
                            numFcn = 0, 
                            numDer = 2)))

classified2 = ddalphaf.classify(c2, test)

print(unlist(classified2))
print(c2$adc.args) 

}
}

\keyword{ functional }
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }
