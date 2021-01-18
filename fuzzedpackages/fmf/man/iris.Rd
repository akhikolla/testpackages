\name{iris}
\alias{iris}
\docType{data}
\title{
Iris Data Set
}
\description{
This is perhaps the best known database to be found in the pattern recognition literature. Fisher's paper is a classic in the field and is referenced frequently to this day. (See Duda & Hart, for example.) 
The data set contains 3 classes of 50 instances each, where each class refers to a type of iris plant. One class is linearly separable from the other 2; the latter are NOT linearly separable from each other. 
}
\usage{data(iris)}
\format{
  A data frame with 150 Instances and 4 attributes (including the class attribute, "Species")
  In this package, the iris dataset has been normalized by the max-min normalization.
}
\details{
Fisher,R.A. "The use of multiple measurements in taxonomic problems" Annual Eugenics, 7, Part II, 179-188 (1936); also in "Contributions to Mathematical Statistics" (John Wiley, NY, 1950). 

Predicted attribute: class of iris plant. 

This is an exceedingly simple domain. 

This data differs from the data presented in Fishers article (identified by Steve Chadwick, spchadwick '@' espeedaz.net ). 
The 35th sample should be: 4.9,3.1,1.5,0.2,"setosa" where the error is in the fourth feature. The 38th sample: 4.9,3.6,1.4,0.1,"setosa" where the errors are in the second and third features. 
}
\source{
Creator: 

R.A. Fisher 

Donor: 

Michael Marshall 
}
\references{
https://archive.ics.uci.edu/ml/datasets/iris
}
\examples{

data(iris)

X = iris[, -1]

y = iris[, 1]
}
\keyword{datasets}
