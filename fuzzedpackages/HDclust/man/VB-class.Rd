\docType{class}
\name{VB-class}
\alias{VB-class}
\alias{VB}
\alias{show,VB-method}

\title{Class "VB" to represent a variable block structure.}
\description{
An S4 class to represent a variable block structure. To create a new
instance of the class, use \code{\link{vb}}.
}
\section{Methods}{
\itemize{
\item \bold{show} signature(object = "VB") : show parameters of variable blocks structure.
\item \bold{getNb} signature(object = "VB") : accessor for 'nb' slot.
\item \bold{getDim} signature(object = "VB") : accessor for 'dim' slot.
\item \bold{getBdim} signature(object = "VB") : accessor for 'bdim' slot.
\item \bold{getNumst} signature(object = "VB") : accessor for 'numst' slot.
\item \bold{getVarorder} signature(object = "VB") : accessor for 'varorder' slot.
}}

\section{Slots}{

\describe{
\item{\code{nb}}{The number of variable blocks.}

\item{\code{dim}}{Dimensionality of the data.}

\item{\code{bdim}}{An integer vector specifying dimensionality of each variable block.}

\item{\code{numst}}{An integer vector specifying the number of mixture models in each 
variable block.}

\item{\code{varorder}}{A list of integer vectors specifying the variable order in 
each variable block.}
}}

