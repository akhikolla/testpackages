\docType{class}
\name{HMMVBBIC-class}
\alias{HMMVBBIC-class}
\alias{HMMVBBIC}
\alias{show,HMMVBBIC-method}
\alias{plot,HMMVBBIC,missing-method}
\alias{getOptHMMVB,HMMVBBIC-method}


\title{Class "HMMVBBIC" to represent results of HMM-VB model selection.}
\description{
An S4 class to represent results of HMM-VB model selection. New instances of the class are created by \code{\link{hmmvbBIC}}.
}
\section{Methods}{
\itemize{
\item \bold{show} signature(object = "HMMVBBIC") : show optimal model.
\item \bold{plot} signature(x = "HMMVBBIC", y = "missing", ...) : plot model selection results (doesn't work for configuration list provided as input to model selection).
\item \bold{getBIC} signature(object = "HMMVBBIC") : accessor for 'BIC' slot.
\item \bold{getLoglikehd} signature(object = "HMMVBBIC") : accessor for 'Loglikehd' slot.
\item \bold{getOptHMMVB} signature(object = "HMMVBBIC") : accessor for 'optHMMVB' slot.
}}
\section{Slots}{

\describe{
\item{\code{BIC}}{A numeric vector specifying calculated BIC values.}
\item{\code{optHMMVB}}{The optimal HMM-VB model with smallest BIC value.}
\item{\code{numst}}{An integer vector specifying the number of mixture components (states) in each variable block for which BIC was calculated. Number of states is the same for all variable blocks.}
}
}
\seealso{
\code{\link{hmmvbBIC}}
}
