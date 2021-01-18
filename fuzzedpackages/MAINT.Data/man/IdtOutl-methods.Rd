\name{IdtOutl-methods}
\docType{methods}
\alias{IdtOutl-methods}
\alias{plot,IdtOutl,missing-method}

\title{Plot method for class IdtOutl in Package `MAINT.Data'}
\description{Plots robust Mahalanobis distances and outlier cut-offs for an object describing potential outliers in a interval-valued data set}

\usage{
  \S4method{plot}{IdtOutl,missing}(x, scale=c("linear","log"), RefDist=getRefDist(x), eta=geteta(x), 
  multiCmpCor=getmultiCmpCor(x), \dots)
}

\arguments{
    \item{x}{An IData object of class IdtOutl describing potential interval-valued ouliters.}

    \item{scale}{The scale of the axis for the robust Mahalanobis distances.}

    \item{RefDist}{The assumed reference distributions used to find cutoffs defining the observations assumed as outliers. Alternatives are \dQuote{ChiSq} and \dQuote{CerioliBetaF} respectivelly for the usual Chi-squared, and the Beta and F distributions proposed by Cerioli (2010). By default uses the one selected in the creation of the object \sQuote{x}.}    

    \item{eta}{Nominal size of the null hypothesis that a given observation is not an outlier. By default uses the one selected in the creation of the object \sQuote{x}.}

    \item{multiCmpCor}{Whether a multicomparison correction of the nominal size (eta) for the outliers tests was performed. Alternatives are: \sQuote{never} -- ignoring the multicomparisons and testing all entities at the \sQuote{eta} nominal level. \sQuote{always} -- testing all n entitites at 1.- (1.-\sQuote{eta}^(1/n)). By default uses the one selected in the creation of the object \sQuote{x}.}   

  \item{\dots}{ Further arguments to be passed to methods.}
}

\references{
Cerioli, A. (2010), Multivariate Outlier Detection with High-Breakdown Estimators.
\emph{Journal of the American Statistical Association} \bold{105} (489), 147--156.\cr 

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr
\emph{Journal of Computational and Graphical Statistics} \bold{14}, 910--927.\cr
}

\keyword{methods}
\keyword{outlier detection}
\keyword{trimmed maximum likelihood estimator}
\keyword{interval data}

\seealso{ \code{\link{getIdtOutl}}, \code{\link{fasttle}}, \code{\link{fulltle}} }



