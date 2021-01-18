\docType{class}
\name{beam-class}
\alias{beam-class}
\alias{print,beam-method}
\alias{show,beam-method}
\alias{summary,beam-method}
\alias{marg,beam-method}
\alias{cond,beam-method}
\alias{mcor,beam-method}
\alias{pcor,beam-method}
\alias{postExpSigma,beam-method}
\alias{postExpOmega,beam-method}
\alias{plotCor,beam-method}
\alias{plotML,beam-method}
\alias{bgraph,beam-method}
\alias{ugraph,beam-method}
\alias{marg}
\alias{cond}
\alias{mcor}
\alias{pcor}
\alias{plotML}
\alias{bgraph}
\alias{ugraph}
\alias{plotCor}
\alias{postExpSigma}
\alias{postExpOmega}

\title{Class beam}

\usage{
  \S4method{print}{beam}(x, ...)

  \S4method{show}{beam}(object)

  \S4method{summary}{beam}(object, ...)

  \S4method{marg}{beam}(object)

  \S4method{cond}{beam}(object)

  \S4method{mcor}{beam}(object)

  \S4method{pcor}{beam}(object)
  
  \S4method{postExpSigma}{beam}(object, vars.method="eb")
  
  \S4method{postExpOmega}{beam}(object, vars.method="eb")

  \S4method{plotML}{beam}(object, ...)

  \S4method{plotCor}{beam}(object, type = object@type, order = 'original', by = "marginal")

  \S4method{bgraph}{beam}(object)

  \S4method{ugraph}{beam}(object)
}

\arguments{
  \item{x}{An object of class \code{beam-class}}
  
  \item{object}{An object of class \code{beam-class}}

  \item{type}{character. Type of correlation to be displayed (marginal, conditional or both)}

  \item{order}{character. Either 'original' or 'clust'. If 'clust' the rows and columns of the correlation matrix are reordered using the cluster memberships obtained by the Louvain clustering algorithm.}

  \item{by}{character. When type ="both" and order = 'clust', specifies whether the clustering has to be performed using the complete weighted marginal or conditional independence graph.}

  \item{vars.method}{character. Method of shrinkage estimation for the variances. Either 'eb', 'mean', 'median' for shrinkage estimation of variance respectively towards an estimated shrinkage target, the mean or the median of the sample variances. Choosing 'none' carries out no shrinkage and uses the sample variances, whereas choosing 'scaled' means that the sample covariance has unit diagonal.}

  \item{...}{further arguments passed to or from other methods.}
  
}

\description{
An S4 class representing the output of the \code{\link{beam}} function.
}

\section{Slots}{
\describe{
  \item{\code{table}}{dat.frame. A data.frame containing marginal and/or partial correlation estimates, Bayes factors and tail probabilities for each edge.}

  \item{\code{deltaOpt}}{numeric. Empirical Bayes estimate of hyperpaprameter delta.}

  \item{\code{alphaOpt}}{numeric. Empirical Bayes estimate of hyperpaprameter alpha.}

  \item{\code{dimX}}{numeric. Dimension of the input data matrix X.}

  \item{\code{type}}{character. Input argument.)}

  \item{\code{varlabs}}{character. Column labels of X.}

  \item{\code{gridAlpha}}{matrix. A matrix containing the log-marginal likelihood of the Gaussian conjugate model as a function of a grid of values of alpha and delta.}

  \item{\code{valOpt}}{numeric. Maximum value of the log-marginal likelihood of the Gaussian conjugate model.}

  \item{\code{return.only}}{character. Input argument.}

  \item{\code{time}}{numeric. Running time (in seconds).}
  
  \item{\code{TinvStdev}}{numeric. Square root of partial variances.}
  
  \item{\code{s}}{numeric. Sample variances.}
}
}

\author{
Gwenael G.R. Leday and Ilaria Speranza
}

