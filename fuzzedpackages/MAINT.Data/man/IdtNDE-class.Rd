\name{IdtNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtNDE-class}

\title{Class IdtNDE}
\description{IdtNDE is a a union of classes \code{\linkS4class{IdtSngNDE}}, \code{\linkS4class{IdtSngNDRE}}, \code{\linkS4class{IdtMxNDE}} and \code{\linkS4class{IdtMxNDRE}}, used for storing the estimation results of Normal modelizations for Interval Data.} 

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\seealso{
  \code{\linkS4class{IdtSngNDE}}, \code{\linkS4class{IdtSngNDRE}}, \code{\linkS4class{IdtMxNDE}}, \code{\linkS4class{IdtMxNDRE}}, \code{\linkS4class{IdtSNDE}}, \code{\linkS4class{IData}}, \code{\link{mle}}, \code{\link{fasttle}}, \code{\link{fulltle}},  \code{\link{MANOVA}}, \code{\link{RobMxtDEst}}
}  

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(coef = "IdtNDE")}: extracts parameter estimates from objects of class IdtNDE}
    \item{stdEr}{\code{signature(x = "IdtNDE")}: extracts standard errors from objects of class IdtNDE}
    \item{vcov}{\code{signature(x = "IdtNDE")}: extracts an estimate of the variance-covariance matrix of the parameters estimators for objects of class IdtNDE}
    \item{mean}{\code{signature(x = "IdtNDE")}: extracts the mean vector estimate from objects of class IdtNDE}
    \item{var}{\code{signature(x = "IdtNDE")}: extracts the variance-covariance matrix estimate from objects of class IdtNDE}
    \item{cor}{\code{signature(x = "IdtNDE")}: extracts the correlation matrix estimate from objects of class IdtNDE}
    \item{sd}{\code{signature(Idt = "IdtNDE")}: extracts the standard deviation estimates from objects of class IdtNDE.}
 }
}

\keyword{classes}
\keyword{Interval Data}



