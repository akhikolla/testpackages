\name{dataf.tecator}
\alias{dataf.tecator}
\alias{tecator}
\docType{data}
\title{
Functional Data Set Spectrometric Data (Tecator)
}
\description{
This dataset is a part of the original one which can be found at
\url{http://stat.cmu.edu}. For each peace of finely chopped meat we observe one spectrometric curve which
corresponds to the absorbance measured at 100 wavelengths. 
The peaces are split according to Ferraty and Vieu (2006) into two classes: with small (<20) and large fat
content obtained by an analytical chemical processing.
}
\usage{
dataf.tecator()
}
\format{
  The functional data as a data structure.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates.
      \describe{
        \item{\code{args}}{\bold{wavelength} - a numeric vector of discretization points from 850 to 1050mm }
        \item{\code{vals}}{\bold{absorbance} - a numeric vector of absorbance values}
       }
    }
    \item{\code{labels}}{The classes of the objects: "small" (<20) and "large" fat content}
  }
}
\author{

Febrero-Bande, M and Oviedo de la Fuente, Manuel
}
\source{

\url{http://stat.cmu.edu}

}
\references{
Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis: theory and practice}. Springer.
}
\seealso{
\code{\link{dataf.*}} for other functional data sets

\code{\link{plot.functional}} for building plots of functional data
}
\examples{
## load the dataset
dataf = dataf.tecator()

## view the classes
unique(dataf$labels)

## access the 5th point of the 2nd object
dataf$dataf[[2]]$args[5]
dataf$dataf[[2]]$vals[5]

## plot the data
\dontrun{
labels = unlist(dataf$labels)
plot(dataf, 
  xlab="Wavelengths", ylab="Absorbances", 
  main=paste("Tecator: < 20 red (", sum(labels == "small"), "),", 
            " >= 20 blue (", sum(labels == "large"), ")", sep=""),
    colors = c("blue", "red"))
}
}
\keyword{datasets}
\keyword{functional}
