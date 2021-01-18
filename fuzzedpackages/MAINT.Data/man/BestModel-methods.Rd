\name{BestModel-methods}
\docType{methods}
\alias{BestModel-methods}
\alias{BestModel,IData-method}   
\alias{BestModel}

\title{Methods for function BestModel in Package `MAINT.Data'}

\description{Selects the best model according to the chosen selection criterion (currently, BIC or AIC) }

\usage{BestModel(ModE,SelCrit=c("IdtCrt","BIC","AIC"))}

\arguments{
  \item{ModE}{An object of class \code{\linkS4class{IdtE}} representing the estimates of a model fitted to a data set of interval-value variables}
  \item{SelCrit}{The model selection criterion. \dQuote{IdtCrt} stands for the criterion originally used in the ModE estimation, while \dQuote{BIC} and \dQuote{AIC} represent respectively the Bayesian and Akaike information criteria.}
}

\value{An integer with the index of the model chosen by the selection criterion}

\keyword{methods}
\keyword{AIC}
\keyword{BIC}
\keyword{interval data}




