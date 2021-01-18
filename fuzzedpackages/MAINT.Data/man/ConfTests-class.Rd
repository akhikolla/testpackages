\name{ConfTests-class}
\Rdversion{1.1}
\docType{class}
\alias{ConfTests-class}
\alias{show,ConfTests-method}

\title{Class "Configuration Tests"}
\description{ConfTests contains a list of the results of statistical likelihood-ratio tests 
that evaluate the goodness-of-fit of restricted models against more general ones. Currently,
the models implemented are those based on the Normal and Skew-Normal distributions, with 
the four alternative variance-covariance matrix configurations.}

\section{Slots}{
  \describe{
    \item{\code{TestRes}:}{List of test results; each element is an object of class LRTest, with the following components: \cr \cr
        \bold{ChiSq}: Value of the Chi-Square statistics corresponding to the performed test. \cr
        \bold{df}: Degrees of freedom of the Chi-Square statistics. \cr
        \bold{pvalue}: p-value of the Chi-Square statistics value, obtained from the Chi-Square distribution with df degrees of freedom. \cr
        \bold{H0logLik}: Logarithm of the Likelihood function under the null hypothesis. \cr
        \bold{H1logLik}: Logarithm of the Likelihood function under the alternative hypothesis. \cr
    }
    \item{\code{RestModels}:}{The restricted model (corresponding to the null hypothesis)}
    \item{\code{FullModels}:}{The full model (corresponding to the alternative hypothesis)}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ConfTests")}: show S4 method for the ConfTests-class}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\seealso{
  \code{\link{mle}}, \code{\linkS4class{IData}}, \code{\linkS4class{LRTest}} 
}

\keyword{classes}
\keyword{interval data}
