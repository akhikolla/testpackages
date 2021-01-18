\docType{class}
\name{HMMVB-class}
\alias{HMMVB-class}
\alias{HMMVB}
\alias{show,HMMVB-method}
\alias{getHmmChain,HMMVB-method}
\alias{getDiagCov,HMMVB-method}
\alias{getVb,HMMVB-method}

\title{Class "HMMVB" to represent a Hidden Markov Model on Variable Blocks (HMM-VB).}
\description{
An S4 class to represent a Hidden Markov Model on Variable Blocks (HMM-VB). 
New instances of the class are created by \code{\link{hmmvbTrain}}.
}
\section{Methods}{
\itemize{
\item \bold{show} signature(object = "HMMVB") : show parameters of the HMM-VB.
\item \bold{getHmmChain} signature(object = "HMMVB") : accessor for 'HmmChain' slot.
\item \bold{getDiagCov} signature(object = "HMMVB") : accessor for 'diagCov' slot.
\item \bold{getBIC} signature(object = "HMMVB") : accessor for 'BIC' slot.
\item \bold{getVb} signature(object = "HMMVB") : accessor for 'VbStructure' slot.
}}
\section{Slots}{

\describe{
\item{\code{VbStructure}}{An object of class 'VB' that contains the variable block structure.}

\item{\code{HmmChain}}{A list of objects of class 'HMM' with trained Hidden Markov Models
for each variable block.}

\item{\code{diagCov}}{A logical value indicating whether or not covariance matrices 
for mixture models are diagonal.}

\item{\code{Loglikehd}}{Loglikelihood value for each data point.}

\item{\code{BIC}}{BIC value for provided variable block structure or optimal BIC value
for found variable block structure.}
}}

