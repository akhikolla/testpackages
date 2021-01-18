\name{dataf.geneexp}
\alias{dataf.geneexp}
\alias{geneexp}
\docType{data}
\title{
Gene Expression Profile Data
}
\description{
A subet of the Drosophila life cycle gene expression data of Arbeitman et al. (2002). The original data set contains 77 gene expression profiles during 58 sequential time points from the embryonic, larval, and pupal periods of the life cycle. The gene expression levels were obtained by a cDNA microarray experiment. 
}
\usage{
dataf.geneexp()
}
\format{
  The functional data as a data structure.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates.
      \describe{
        \item{\code{args}}{\bold{Time} - a numeric vector of time periods}
        \item{\code{vals}}{\bold{Gene Expression Level} - a numeric vector}
       }
    }
    \item{\code{labels}}{ 
    Biological classifications identified in Arbeitman et al.(2002) (1 = transient early zygotic genes; 2 = muscle-specific genes; 3 = eye-specific genes. )}
  }
}
\source{
Chiou, J.-M. and Li, P.-L. Functional clustering and identifying substructures of longitudinal data, J. R. Statist. Soc. B, Volume 69 (2007), 679-699

Arbeitman, M.N., Furlong, E.E.M., Imam,F., Johnson, E., Null,B.H., Baker,B.S., Krasnow, M.A., Scott,M.P., Davis,R.W. and White,K.P. (2002) Gene expression during the life cycle of Drosophila melanogaster. Science, 297, 2270-2274.

}
\seealso{
\code{\link{dataf.*}} for other functional data sets

\code{\link{plot.functional}} for building plots of functional data
}
\examples{
## load the dataset
dataf = dataf.geneexp()

## view the classes
unique(dataf$labels)

## access the 5th point of the 2nd object
dataf$dataf[[2]]$args[5]
dataf$dataf[[2]]$vals[5]

## plot the data
\dontrun{
labels = unlist(dataf$labels)
plot(dataf, 
  xlab="Time", ylab="Gene Expression Level", 
  main=paste0("Gene Expression:  1 red (", sum(labels == 1), "), ", 
            "2 green (", sum(labels == 2), "), ", 
            "3 blue (", sum(labels == 3), ")"),
  colors = c("red", "green", "blue"))
}
}
\keyword{datasets}
\keyword{functional}
