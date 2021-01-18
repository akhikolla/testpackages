\name{LRTest-class}
\Rdversion{1.1}
\docType{class}
\alias{LRTest-class}
\alias{show,LRTest-method}

\title{Class LRTest}
\description{LRTest contains the results of likelihood ratio tests
}

\section{Slots}{
  \describe{
   \item{\code{ChiSq}:}{Value of the Chi-Square statistics corresponding to the performed test}
    \item{\code{df}:}{Degrees of freedom of the Chi-Square statistics}
    \item{\code{pvalue}:}{p-value of the Chi-Square statistics value, obtained from the Chi-Square distribution with df degrees of freedom}
    \item{\code{H0logLik}:}{Logarithm of the Likelihood function under the null hypothesis}
    \item{\code{H1logLik}:}{Logarithm of the Likelihood function under the alternative hypothesis}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "LRTest")}: show S4 method for the LRTest-class
 }
	 }
}
\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\seealso{
 \code{\link{mle}}, \code{\linkS4class{IData}}, 
\code{\linkS4class{ConfTests}}, \code{\link{MANOVA}}
}

\keyword{classes}
\keyword{interval data}
\keyword{likelihood ratio test}
