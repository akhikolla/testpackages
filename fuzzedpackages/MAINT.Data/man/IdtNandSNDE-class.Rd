\name{IdtNandSNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtNandSNDE-class}

\title{Class IdtNandSNDE}
\description{IdtNandSNDE is a union of classes \code{\linkS4class{IdtSngNandSNDE}} and \code{\linkS4class{IdtMxNandSNDE}}, used for storing the estimation results of Normal and Skew-Normal modelisations for Interval Data.} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(coef = "IdtNandSNDE")}: extracts parameter estimates from objects of class IdtNandSNDE}
    \item{stdEr}{\code{signature(x = "IdtNandSNDE")}: extracts standard errors from objects of class IdtNandSNDE}
    \item{vcov}{\code{signature(x = "IdtNandSNDE")}: extracts an estimate of the variance-covariance matrix of the parameters estimators for objects of class IdtNandSNDE}
    \item{mean}{\code{signature(x = "IdtNandSNDE")}: extracts the mean vector estimate from objects of class IdtNandSNDE}
    \item{var}{\code{signature(x = "IdtNandSNDE")}: extracts the variance-covariance matrix estimate from objects of class IdtNandSNDE}
    \item{cor}{\code{signature(x = "IdtNandSNDE")}: extracts the correlation matrix estimate from objects of class IdtNandSNDE}
 }
}

\seealso{
  \code{\linkS4class{IData}}, \code{\link{mle}}, \code{\link{fasttle}}, \code{\link{fulltle}}, 
  \code{\link{MANOVA}}, \code{\link{RobMxtDEst}}, \code{\linkS4class{IdtSngNandSNDE}}, \code{\linkS4class{IdtMxNandSNDE}}  
}

\keyword{classes}
\keyword{Interval Data}



