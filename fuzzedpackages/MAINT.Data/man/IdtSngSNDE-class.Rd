\name{IdtSngSNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtSngSNDE-class}

\title{Class IdtSngSNDE}
\description{Contains the results of a single class maximum likelihood estimation for the Skew-Normal distribution, with the four different possible variance-covariance configurations.
}

\section{Slots}{
  \describe{
   \item{\code{CovConfCases}:}{List of the considered configurations}
    \item{\code{ModelNames}:}{The model acronym, indicating the model type (currently, N for Normal and SN for Skew-Normal), 
and the configuration Case (C1 to C4) for the covariance matrix}
    \item{\code{ModelNames}:}{Inherited from class \code{\linkS4class{IdtE}}. The model acronym formed by a "SN", indicating a skew-Normal model, followed by the configuration (Case 1 through Case 4)}
    \item{\code{ModelType}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the model; always set to "SkewNormal" in objects of the IdtSngSNDE class}    
\item{\code{ModelConfig}:}{Inherited from class \code{\linkS4class{IdtE}}. Configuration case of the variance-covariance matrix: Case 1 through Case 4}    
    \item{\code{NIVar}:}{Inherited from class \code{\linkS4class{IdtE}}. Number of interval variables}
    \item{\code{SelCrit}:}{Inherited from class \code{\linkS4class{IdtE}}. The model selection criterion; currently, AIC and BIC are implemented}
    \item{\code{logLiks}:}{Inherited from class \code{\linkS4class{IdtE}}. The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the AIC criterion }
    \item{\code{BICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the BIC criterion }
    \item{\code{BestModel}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the best model according to the chosen selection criterion}
    \item{\code{SngD}:}{Inherited from class \code{\linkS4class{IdtE}}. Boolean flag indicating whether a single or a mixture of distribution were estimated. Always set to TRUE in objects of class IdtSngSNDE }
  }
}

\section{Extends}{
Class \code{\linkS4class{IdtSngDE}}, directly.
Class \code{\linkS4class{IdtE}}, by class \code{\linkS4class{IdtSngDE}}, distance 2.
}

\section{Methods}{
No methods defined with class IdtSngSNDE in the signature.
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Azzalini, A. and Dalla Valle, A. (1996), The multivariate skew-normal distribution. \emph{Biometrika} \bold{83}(4), 715--726.\cr

Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\seealso{
   \code{\link{mle}}, \code{\linkS4class{IData}}, \code{\linkS4class{IdtSngNDE}}, \code{\linkS4class{IdtMxSNDE}}  
}

\keyword{classes}
\keyword{interval data}
