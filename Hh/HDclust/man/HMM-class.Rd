\docType{class}
\name{HMM-class}
\alias{HMM-class}
\alias{HMM}
\alias{show,HMM-method}
\alias{getPrenumst,HMM-method}
\alias{getHmmParam,HMM-method}
\title{Class "HMM" to represent parameters associated with a variable block in the HMM-VB}
\description{
An S4 class to represent the model parameters associated with one variable block in the HMM-VB.
For brevity, we call this part of HMM-VB, specific to a particular variable block, an "HMM" for the block. New instances of the class are created by \code{\link{hmmvbTrain}}.
}
\section{Methods}{
\itemize{
\item \bold{show} signature(object = "HMM") : show parameters of the HMM object.
\item \bold{getPrenumst} signature(object = "HMM") : accessor for 'prenumst' slot.
\item \bold{getHmmParam} signature(object = "HMM") : accessor for parameters of the HMM object. This function outputs a list with means, covariance matrices, inverse covarince matrices and logarithms of the determinants of the covariance matrices for all states of the HMM.
}}
\section{Slots}{

\describe{
\item{\code{dim}}{Dimensionality of the data in HMM.}

\item{\code{numst}}{An integer vector specifying the number of HMM states.}

\item{\code{prenumst}}{An integer vector specifying the number of states of previous 
variable block HMM.}

\item{\code{a00}}{Probabilities of HMM states.}

\item{\code{a}}{Transition probability matrix from states in the previous variable block
to the states in the current one.}

\item{\code{mean}}{A numerical matrix with state means. \emph{k}th row corresponds to the
\emph{k}th state.}

\item{\code{sigma}}{A list containing the covariance matrices of states.}

\item{\code{sigmaInv}}{A list containing the inverse covariance matrices of states.}

\item{\code{sigmaDetLog}}{A vector with \eqn{log(|sigma|)} for each state.}
}}

