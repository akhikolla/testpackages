\docType{class}
\name{HMMVBclust-class}
\alias{HMMVBclust-class}
\alias{HMMVBclust}
\alias{show,HMMVBclust-method}
\alias{plot,HMMVBclust,missing-method}
\alias{getClustParam,HMMVBclust-method}
\alias{getClsid,HMMVBclust-method}
\alias{getSize,HMMVBclust-method}


\title{Class "HMMVBclust" to represent clustering results with Hidden Markov Model on variable block structure.}
\description{
An S4 class to represent a clustering result based on HMM-VB. New instances of the class are created by \code{\link{hmmvbClust}}.
}
\section{Methods}{
\itemize{
\item \bold{show} signature(object = "HMMVBclust") : show clustering results based on HMM-VB.
\item \bold{plot} signature(x = "HMMVBclust", y = "missing", method = "t-sne", ...) : plot clustering results. 'method' controls the visualization algorithm. Two algorithms are supported:  method = 'PCA' plots the data using 2 component PCA space; and method = 't-SNE' plots the data using 2 component t-SNE space. Default setting is t-SNE.
\item \bold{getClustParam} signature(object = "HMMVBclust") : accessor for 'clustParam' slot.
\item \bold{getLoglikehd} signature(object = "HMMVBclust") : accessor for 'Loglikehd' slot.
\item \bold{getClsid} signature(object = "HMMVBclust") : accessor for 'clsid' slot.
\item \bold{getSize} signature(object = "HMMVBclust") : accessor for 'size' slot.
}}
\section{Slots}{

\describe{
\item{\code{data}}{The input data matrix}

\item{\code{clustParam}}{A list with cluster parameters:

\describe{
\item{ncls}{The number of clusters (same as the number of modes)}
\item{mode}{A numeric matrix with cluster modes. \emph{k}th row of the matrix 
stores coordinates of the \emph{k}th mode.}
\item{ndseq}{The number of distinct Viterbi sequences for the dataset}
\item{vseqid}{An integer vector representing the map between Viterbi 
sequences and clusters. \emph{k}th value in the vector stores cluster id for
\emph{k}th Viterbi sequence.}
\item{vseq}{A list with integer vectors representing distinct Viterbi sequences
for the dataset}
\item{sigma}{A numeric vector with the dataset variance}
}}

\item{\code{clsid}}{An integer vector with cluster ids.}
\item{\code{Loglikehd}}{Loglikelihood value for each data point.}
\item{\code{size}}{An integer vector with cluster sizes.}
}}

