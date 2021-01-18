\docType{class}
\name{beam.select-class}
\alias{beam.select-class}
\alias{print,beam.select-method}
\alias{show,beam.select-method}
\alias{summary,beam.select-method}
\alias{marg,beam.select-method}
\alias{cond,beam.select-method}
\alias{mcor,beam.select-method}
\alias{pcor,beam.select-method}
\alias{plotAdj,beam.select-method}
\alias{plotML,beam.select-method}
\alias{bgraph,beam.select-method}
\alias{ugraph,beam.select-method}
\alias{plotAdj}

\title{Class beam.select}

\usage{
  \S4method{print}{beam.select}(x, ...)

  \S4method{show}{beam.select}(object)

  \S4method{summary}{beam.select}(object, ...)
  
  \S4method{marg}{beam.select}(object)

  \S4method{cond}{beam.select}(object)

  \S4method{mcor}{beam.select}(object)
  
  \S4method{pcor}{beam.select}(object)

  \S4method{plotML}{beam.select}(object, ...)

  \S4method{plotAdj}{beam.select}(object, type=object@type, order = "original")

  \S4method{bgraph}{beam.select}(object)

  \S4method{ugraph}{beam.select}(object)

}

\arguments{
  \item{x}{An object of class \code{beam.select-class}}
  
  \item{object}{An object of class \code{beam.select-class}}

  \item{type}{character. Type of correlation to be displayed (marginal, conditional or both)}

  \item{order}{character. Either 'original' or 'clust'. If 'clust' the rows and columns of the adjacency matrix are reordered using the cluster memberships obtained by the Louvain clustering algorithm.}

  \item{...}{further arguments passed to or from other methods.}
  
}

\description{
An S4 class representing the output of the \code{\link{beam.select}} function.
}

\section{Slots}{
\describe{
  \item{\code{marginal}}{data.frame. A data.frame containing the marginal correlation estimates, Bayes factors and tail probabilities for the selected edges only.}

  \item{\code{conditional}}{data.frame. A data.frame containing the partial correlation estimates, Bayes factors and tail probabilities for the selected edges only.}

  \item{\code{dimX}}{numeric. Dimension of the imput data matrix X.}

  \item{\code{type}}{character. Input type (marginal, conditional or both)}

  \item{\code{varlabs}}{character. Column labels of X.}

  \item{\code{alphaOpt}}{numeric. Empirical Bayes estimates of hyperpaprameter alpha.}

  \item{\code{gridAlpha}}{matrix. A matrix containing the log-marginal likelihood of the Gaussian conjugate model as a function of a grid of values of alpha and delta.}

  \item{\code{valOpt}}{numeric. Maximum value of the log-marginal likelihood of the Gaussian conjugate model}

  \item{\code{method}}{character. Input method.}

  \item{\code{thres}}{numeric. Input threshold}
}
}

\author{
Gwenael G.R. Leday and Ilaria Speranza
}

