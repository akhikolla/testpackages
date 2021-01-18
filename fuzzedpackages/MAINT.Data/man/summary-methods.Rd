\name{summary-methods}
\docType{methods}
\alias{summary,IdtMclust-method}
\title{ IdtMclust summary method}

\description{summary methods for the classe IdtMclust defined in Package `MAINT.Data'.}

\usage{
 \S4method{summary}{IdtMclust}(object, parameters = FALSE, classification = FALSE, model = "BestModel", \dots)
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{IdtMclust}} representing the results of fitting Gaussian mixtures to interval data objects}
  \item{parameters}{A boolean flag indicating if the parameter estimates of the optimal mixture should be displayed}
  \item{classification}{A boolean flag indicating if the crisp classification resulting from the optimal mixture should be displayed}
  \item{model}{A character vector specifying the the model whose solution is to be displayed.}
  \item{\dots}{Other named arguments.}
}

\seealso{
  \code{\link{Idtmclust}}, \code{\linkS4class{IdtMclust}}, \code{\link{plotInfCrt}},  \code{\link{pcoordplot}} 
}


\keyword{methods}
\keyword{interval data}




