\name{beam.select}
\alias{beam.select}

\title{Edge selection with multiple testing and error control}

\usage{
beam.select(object, thres = 0.1, method = "BH",
return.only = c(object@return.only, "adj"))
}

\arguments{
  \item{object}{An object of class \code{\link{beam-class}}.}

  \item{thres}{numeric. Threshold to be applied on adjusted tail probabilities.}

  \item{method}{character. Method to use for multiple comparison adjustment of tail probabilities.}

  \item{return.only}{character. Quantities to be returned.}
}

\value{
An object of class \code{\link{beam.select-class}}
}

\description{
Infer graphical structures by multiple testing
}

\details{
The argument \code{method} allows to adjust the tail probabilities obtained from the null distributions of
the Bayes factors for multiple comparisons. Possible choices are: "holm", "bonferroni", "BH", "BY" and "HC".
Apart from "HC", these are passed onto the R function \code{p.adjust}
from package \pkg{stats} and we refer the user to its documentation for details. The method "HC" provides an
optimal decision threshold based on the Higher Criticism score which is computed using the R function \code{hc.thresh}
from package \pkg{fdrtool}. Again, we refer to the associated documentation for details.

The argument \code{return.only} allows to decide which quantities have to be in the output: it could be any subvector of c('cor', 'BF', 'prob', 'adj') (provided that the requested quantities have been computed in the beam object, except for adjusted probabilities). It can also be set to NULL: in this case, only the selected edges will be returned without any additional information. The default value for this argument are the columns present in the beam object plus the adjusted probabilities.
}

\author{
Gwenael G.R. Leday and Ilaria Speranza
}

\references{
Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 430-449.\cr
Goeman, J. J., & Solari, A. (2014). Multiple hypothesis testing in genomics. Statistics in medicine, 33(11), 1946-1978.\cr
Donoho, D., & Jin, J. (2015). Higher criticism for large-scale inference, especially for rare and weak effects. Statistical Science, 30(1), 1-25.\cr
Klaus, B., & Strimmer, K. (2012). Signal identification for rare and weak features: higher criticism or false discovery rates?. Biostatistics, 14(1), 129-143.
}

