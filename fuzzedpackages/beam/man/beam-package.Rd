\name{beam-package}
\alias{beam-package}
\docType{package}

\title{Fast Bayesian Inference in Large Gaussian Graphical Models}

\description{
The package enables inference of marginal and conditional dependencies from high-dimensional using the method of Leday and Richardson (2019). Inference is carried out by multiple testing of hypotheses about pairwise (marginal or conditional) independence using closed-form Bayes factors. Exact tail probabilities are obtained from the null distributions of the Bayes factors to help address the multiplicity problem and control desired error rates for incorrect edge inclusion. The method is computationally very efficient and allows to address problems with hundreds or thousands of variables.
}

\details{
\strong{1.} The main function of the package is \code{\link{beam}} which carries out shrinkage estimation of the (inverse) covariance and compute the (scaled) Bayes factors as well as the tail probabilities (p-values). The function returns an (S4) object of class \code{\link{beam-class}} that is associated with the following methods:

- \code{\link{summary,beam-method}}:\cr
provides a summary of inferred (marginal and/or conditional) associations.

- \code{\link{marg,beam-method}}:\cr
returns a data.frame with marginal correlations, Bayes factors and/or tail probabilities.

- \code{\link{cond,beam-method}}:\cr
returns a data.frame with partial correlations, Bayes factors and/or tail probabilities.

- \code{\link{mcor,beam-method}}:\cr
return marginal correlation matrix (scaled posterior expectation of the covariance matrix).

- \code{\link{pcor,beam-method}}:\cr
return partial correlation matrix (scaled posterior expectation of the inverse covariance matrix).

- \code{\link{plotML,beam-method}}:\cr
plot log-marginal likelihood of the Gaussian conjugate model as a function of shrinkage parameter.

- \code{\link{plotCor,beam-method}}:\cr
plot heatmap of marginal (upper triangle) and/or partial (lower triangle) correlation estimates.\cr

\strong{2.} The function \code{\link{beam.select}} takes as input an object of class \code{\link{beam-class}} and carries out edge selection by multiple testing of hypotheses about pairwise (marginal or conditional) independence. The function helps address the multiplicity problem and control different types of error rates (e.g. false discovery rate, family-wise error rate, ...). \code{\link{beam.select}} returns an (S4) object of class \code{\link{beam.select-class}} that is associated with the following methods:

- \code{\link{summary,beam.select-method}}:\cr
provides a summary of inferred (marginal or conditional) associations.

- \code{\link{marg,beam.select-method}}:\cr
returns a data.frame with marginal correlations, Bayes factors and/or tail probabilities for selected edges.

- \code{\link{cond,beam.select-method}}:\cr
returns a data.frame with partial correlations, Bayes factors and/or tail probabilities for selected edges.

- \code{\link{bgraph,beam.select-method}}:\cr
return an \code{igraph} object containing the marginal (in)dependence graph.

- \code{\link{ugraph,beam.select-method}}:\cr
return an \code{igraph} object containing the conditional (in)dependence graph.
}

\author{
Authors: Gwenael G.R. Leday and Ilaria Speranza

Maintainer: Gwenael G.R. Leday <gwenael.leday@mrc-bsu.cam.ac.uk>
}

\references{
Leday, G.G.R. and Richardson, S. (2019). Fast Bayesian inference in large Gaussian graphical models. \emph{Biometrics}. 75(4), 1288--1298.
}
