\name{IdtMclust-class}
\Rdversion{1.1}
\docType{class}

\alias{IdtMclust-class}

\alias{parameters}
\alias{pro}
\alias{classification}
\alias{SelCrit}
\alias{Hmcdt}
\alias{BestG}
\alias{BestC}
\alias{PostProb}

\alias{parameters,IdtMclust-method}
\alias{pro,IdtMclust-method}
\alias{mean,IdtMclust-method}
\alias{var,IdtMclust-method}
\alias{cor,IdtMclust-method}
\alias{classification,IdtMclust-method}
\alias{SelCrit,IdtMclust-method}
\alias{Hmcdt,IdtMclust-method}
\alias{BestG,IdtMclust-method}
\alias{BestC,IdtMclust-method}
\alias{PostProb,IdtMclust-method}
\alias{show,IdtMclust-method}
\alias{BIC,IdtMclust-method}
\alias{AIC,IdtMclust-method}
\alias{logLik,IdtMclust-method}

\title{Class IdtMclust}
\description{IdtMclust contains the results of fitting mixtures of Gaussian distributions to interval data represented by objects of class \code{\linkS4class{IData}}.}

\section{Slots}{
  \describe{
    \item{\code{call}:}{The matched call that created the IdtMclust object} 
    \item{\code{data}:}{The IData data object} 
    \item{\code{NObs}:}{Number of entities under analysis (cases)}
    \item{\code{NIVar}:}{Number of interval variables}
    \item{\code{SelCrit}:}{The model selection criterion; currently, AIC and BIC are implemented}
    \item{\code{Hmcdt}:}{Indicates whether the optimal model corresponds to a homoscedastic (TRUE) or a hetereocedasic (FALSE) setup}
    \item{\code{BestG:}}{The optimal number of mixture components.}
    \item{\code{BestC:}}{The configuration case of the variance-covariance matrix in the optimal model}
    \item{\code{logLiks}:}{The logarithms of the likelihood function for the different models tried}
    \item{\code{logLik}:}{The logarithm of the likelihood function for the optimal model}
    \item{\code{AICs}:}{The values of the AIC criterion for the different models tried}
    \item{\code{aic}:}{The value of the AIC criterion for the he optimal model}
    \item{\code{BICs}:}{The values of the BIC criterion for the different models tried}
    \item{\code{bic}:}{The value of the BIC criterion for the he optimal model}
    \item{\code{parameters}}{
      A list with the following components:
      \describe{
        \item{pro}{A vector whose \emph{k}th component is the mixing proportion for the \emph{k}th component of the mixture model.}
        \item{mean}{The mean for each component. If there is more than one component, this is a matrix whose kth column is the mean of the \emph{k}th component of the mixture model.}
        \item{covariance}{A three-dimensional array with the covariance estimates. If Hmcdt is FALSE (heteroscedastic setups) the third dimension levels run through the BestG mixture components, with one different covariance matrix for each level. Otherwise (homoscedastic setups), there is only one covariance matrix and the size of the third dimension equals one.}
      }
    }
    \item{\code{z:}}{A matrix whose \emph{[i,k]}th entry is the probability that observation \emph{i} in the test data belongs to the \emph{k}th class.}
    \item{\code{classification:}}{The classification corresponding to \code{z}, i.e. \code{map(z)}.}
    \item{\code{allres:}}{A list with the detailed results for all models fitted.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "IdtMclust")}: show S4 method for the IdtMclust-class}
    \item{summary}{\code{signature(object = "IdtMclust")}: summary S4 method for the IdtMclust-class}
    \item{parameters}{\code{signature(x = "IdtMclust")}: retrieves the value of the parameter estimates for the obtained partition}
    \item{pro}{\code{signature(x = "IdtMclust")}: retrieves the value of the estimated mixing proportions for the obtained partition}
    \item{mean}{\code{signature(x = "IdtMclust")}: retrieves the value of the component means for the obtained partition}
    \item{var}{\code{signature(x = "IdtMclust")}: retrieves the value of the estimated covariance matrices for the obtained partition}
    \item{cor}{\code{signature(x = "IdtMclust")}: retrieves the value of the estimated correlation matrices}
    \item{classification}{\code{signature(x = "IdtMclust")}: retrieves the individual class assignments for the obtained partition}
    \item{SelCrit}{\code{signature(x = "IdtMclust")}: retrieves a string specifying the criterion used to find the best model and partition}
    \item{Hmcdt}{\code{signature(x = "IdtMclust")}: returns TRUE if an homecedastic model has been assumed, and FALSE otherwise}
    \item{BestG}{\code{signature(x = "IdtMclust")}: returns the number of components selectd}
    \item{BestC}{\code{signature(x = "IdtMclust")}: retruns the covariance configuration selected}
    \item{PostProb}{\code{signature(x = "IdtMclust")}: retrieves the estimates of the individual posterir probabilities for the obtained partition}
    \item{BIC}{\code{signature(x = "IdtMclust")}: returns the value of the BIC criterion}
    \item{AIC}{\code{signature(x = "IdtMclust")}: returns the value of the AIC criterion}
    \item{logLik}{\code{signature(x = "IdtMclust")}: returns the value of the log-likelihood}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Brito, P., Duarte Silva, A. P. and Dias, J. G. (2015), Probabilistic Clustering of Interval Data. \emph{Intelligent Data Analysis} \bold{19}(2), 293--313.
}

\seealso{
  \code{\link{Idtmclust}}, \code{\link{plotInfCrt}},  \code{\link{pcoordplot}} 
}

\keyword{classes}
\keyword{interval data}
\keyword{model based clustering}

