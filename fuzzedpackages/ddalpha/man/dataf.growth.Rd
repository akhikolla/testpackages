\name{dataf.growth}
\alias{dataf.growth}
\alias{growth}
\docType{data}
\title{
Berkeley Growth Study Data
}
\description{
The data set contains the heights of 39 boys and 54 girls from age 1 to 18 and the ages at which they were collected.
}
\usage{
dataf.growth()
}
\format{
  The functional data as a data structure.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates 
      \describe{
        \item{\code{args}}{\bold{age} - a numeric vector of length 31 giving the ages at which the heights were measured}
        \item{\code{vals}}{\bold{height} - a numeric vector of heights in centimeters of 39 boys and 54 girls}
       }
    }
    \item{\code{labels}}{The classes of the objects: boy, girl}
  }
}
\details{
The ages are not equally spaced.
}
\source{
Ramsay, James O., and Silverman, Bernard W. (2006), Functional Data Analysis, 2nd ed., Springer, New York.

Ramsay, James O., and Silverman, Bernard W. (2002), Applied Functional Data Analysis, Springer, New York, ch. 6.

Tuddenham, R. D., and Snyder, M. M. (1954) "Physical growth of California boys and girls from birth to age 18", University of California Publications in Child Development, 1, 183-364.

}
\seealso{
\code{\link{dataf.*}} for other functional data sets

\code{\link{plot.functional}} for building plots of functional data
}
\examples{
## load the Growth dataset
dataf = dataf.growth()

## view the classes
unique(dataf$labels)

## access the 5th point of the 2nd object
dataf$dataf[[2]]$args[5]
dataf$dataf[[2]]$vals[5]

## plot the data
\dontrun{
  labels = unlist(dataf$labels)
  plot(dataf, 
    main = paste("Growth: girls red (", sum(labels == "girl"), "),", 
                      " boys blue (", sum(labels == "boy"), ")", sep=""),
    xlab="Year", ylab="Height, cm",
    colors = c("blue", "red")   # in alphabetical order of class labels   
    )
}
}
\keyword{datasets}
\keyword{functional}
