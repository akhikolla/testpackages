\name{IdtMxE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMxE-class}

\title{Class IdtMxE}
\description{IdtMxE extends the \code{\linkS4class{IdtE}} class, assuming that the data can be characterized by  a mixture of distributions, for instances considering partitions of entities into different groups.
}

\section{Slots}{
  \describe{
    \item{\code{grouping}:}{Factor indicating the group to which each observation belongs to}
    \item{\code{ModelNames}:}{Inherited from class \code{\linkS4class{IdtE}}. The model acronym, indicating the model type (currently, N for Normal and SN for Skew-Normal), and the configuration (Case 1 through Case 4)}
    \item{\code{ModelType}:}{Inherited from class \code{\linkS4class{IdtE}}. Indicates the model; currently, Gaussian or Skew-Normal distributions are implemented.}
    \item{\code{ModelConfig}:}{Inherited from class \code{\linkS4class{IdtE}}. Configuration of the variance-covariance matrix: Case 1 through Case 4}
    \item{\code{NIVar}:}{Inherited from class \code{\linkS4class{IdtE}}. Number of interval variables}
    \item{\code{SelCrit}:}{Inherited from class \code{\linkS4class{IdtE}}. The model selection criterion; currently, AIC and BIC are implemented}
    \item{\code{logLiks}:}{Inherited from class \code{\linkS4class{IdtE}}. The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the AIC criterion }
    \item{\code{BICs}:}{Inherited from class \code{\linkS4class{IdtE}}. Value of the BIC criterion }
    \item{\code{BestModel}:}{Inherited from class \code{\linkS4class{IdtE}}. Bestmodel indicates the best model according to the chosen selection criterion}
    \item{\code{SngD}:}{Inherited from class \code{\linkS4class{IdtE}}. Boolean flag indicating whether a single or a mixture of distribution were estimated. Always set to FALSE in objects of class "IdtMxE" }
    \item{\code{Ngrps}:}{Number of mixture components}
  }
}

\section{Extends}{
Class \code{\linkS4class{IdtE}}, directly.
}

\section{Methods}{
No methods defined with class "IdtMxE" in the signature.
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\seealso{
  \code{\linkS4class{IdtE}}, \code{\linkS4class{IdtSngDE}}, \code{\linkS4class{IData}}, \code{\link{MANOVA}}, \code{\link{RobMxtDEst}}  
}

\keyword{classes}
\keyword{interval data}

