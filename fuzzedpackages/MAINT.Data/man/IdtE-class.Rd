\name{IdtE-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtE-class}
\alias{IdtSngDE-class}
\alias{summary,IdtE-method}
\alias{show,IdtE-method}
\alias{BestModel,IdtE-method}
\alias{testMod,IdtE-method}
\alias{sd,IdtE-method}
\alias{AIC,IdtE-method}
\alias{BIC,IdtE-method}
\alias{logLik,IdtE-method}

\title{Class IdtE}
\description{IdtE contains estimation results for the models assumed for single distributions, or mixtures of distributions, underlying data sets of interval-valued entities.}

\section{Slots}{
  \describe{
    \item{\code{ModelNames}:}{The model acronym, indicating the model type (currently, N for Normal and SN for Skew-Normal), 
and the configuration (Case 1 through Case 4)}
    \item{\code{ModelType}:}{Indicates the model; currently, Gaussian or Skew-Normal distributions are implemented}
    \item{\code{ModelConfig}:}{Configuration of the variance-covariance matrix: Case 1 through Case 4}
    \item{\code{NIVar}:}{Number of interval variables}
    \item{\code{SelCrit}:}{The model selection criterion; currently,  AIC and BIC are implemented}
    \item{\code{logLiks}:}{The logarithms of the likelihood function for the different cases}
    \item{\code{AICs}:}{Value of the AIC criterion }
    \item{\code{BICs}:}{Value of the BIC criterion }
    \item{\code{BestModel}:}{Bestmodel indicates the best model according to the chosen selection criterion }
    \item{\code{SngD}:}{Boolean flag indicating whether a single or a mixture of distribution were estimated }
  }
}

\section{Methods}{
  \describe{
    \item{BestModel}{\code{signature(Idt = "IdtE")}: Selects the best model according to the chosen selection criterion (currently, AIC or BIC) }
    \item{show}{\code{signature(object = "IdtE")}: show S4 method for the IDtE-class}
    \item{summary}{\code{signature(object = "IdtE")}: summary S4 method for the IDtE-class}
    \item{testMod}{\code{signature(Idt = "IdtE")}: Performs statistical likelihood-ratio tests that evaluate the goodness-of-fit of a nested model against a more general one.}
    \item{sd}{\code{signature(Idt = "IdtE")}: extracts the standard deviation estimates from objects of class IdtE.}
    \item{AIC}{\code{signature(Idt = "IdtE")}: extracts the value of the Akaike Information Criterion from objects of class IdtE.}
    \item{BIC}{\code{signature(Idt = "IdtE")}: extracts the value of the Bayesian Information Criterion from objects of class IdtE.}
    \item{logLik}{\code{signature(Idt = "IdtE")}: extracts the value of the maximised log-likelihood from objects of class IdtE.}
    \item{mean}{\code{signature(x = "IdtE")}: extracts the mean vector estimate from objects of class IdtE}
    \item{var}{\code{signature(x = "IdtE")}: extracts the variance-covariance matrix estimate from objects of class IdtE}
    \item{cor}{\code{signature(x = "IdtE")}: extracts the correlation matrix estimate from objects of class IdtE}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.
}

\seealso{
  \code{\link{mle}}, \code{\link{fasttle}}, \code{\link{fulltle}}, \code{\link{MANOVA}}, \code{\link{RobMxtDEst}}, 
  \code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}

