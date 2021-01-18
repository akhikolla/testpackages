\name{IdtMxSNDE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMxSNDE-class}

\title{Class IdtMxSNDE}
\description{IdtMxSNDE contains the results of a mixture model estimation for the Skew-Normal model, with the four different possible variance-covariance configurations.
}

\section{Slots}{
  \describe{
     \item{\code{Hmcdt}:}{Indicates whether we consider an homoscedastic location model (TRUE) or a general model (FALSE)}
     \item{\code{CovConfCases}:}{List of the considered configurations}
     \item{\code{grouping}:}{Inherited from class \code{\linkS4class{IdtMxE}}. Factor indicating the group to which each observation belongs to}
     \item{\code{ModelNames}:}{Inherited from class \code{\linkS4class{IdtE}}. The model acronym, indicating the model type (currently, N for Normal and SN for Skew-Normal), and the configuration (Case 1 through Case 4)}
     \item{\code{ModelType}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the model; currently, Gaussian or Skew-Normal distributions are implemented}    
     \item{\code{ModelConfig}:}{Inherited from class \code{\linkS4class{IdtE}}. Configuration case of the variance-covariance matrix: Case 1 through Case 4}    
     \item{\code{NIVar}:}{Inherited from class \code{\linkS4class{IdtE}}. Number of interval variables}
     \item{\code{SelCrit}:}{Inherited from class \code{\linkS4class{IdtE}}. The model selection criterion; currently, AIC and BIC are implemented}
     \item{\code{logLiks}:}{Inherited from class \code{\linkS4class{IdtE}}. The logarithms of the likelihood function for the different cases}
     \item{\code{AICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the AIC criterion }
     \item{\code{BICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the BIC criterion }
     \item{\code{BestModel}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the best model according to the chosen selection criterion}
    \item{\code{SngD}:}{Inherited from class \code{\linkS4class{IdtE}}. Boolean flag indicating whether a single or a mixture of distribution were estimated. Always set to FALSE in objects of class \code{\linkS4class{IdtMxSNDE}} }
    \item{\code{Ngrps}:}{Inherited from class \code{\linkS4class{IdtMxE}}. Number of mixture components}
  }  
}

\section{Extends}{
Class \code{\linkS4class{IdtMxE}}, directly.
Class \code{\linkS4class{IdtE}}, by class \code{\linkS4class{IdtMxE}}, distance 2.
}

\section{Methods}{
No methods defined with class IdtMxSNDE in the signature.
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Azzalini, A. and Dalla Valle, A. (1996), The multivariate skew-normal distribution. \emph{Biometrika} \bold{83}(4), 715--726.\cr

Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\seealso{
  \code{\linkS4class{IdtE}}, \code{\linkS4class{IdtMxE}}, \code{\linkS4class{IdtSngSNDE}}, \code{\link{MANOVA}}, \code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}

