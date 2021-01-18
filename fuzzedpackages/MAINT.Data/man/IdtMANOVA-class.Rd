\name{IdtMANOVA-class}
\Rdversion{1.1}
\docType{class}
\alias{IdtMANOVA-class}
\alias{IdtClMANOVA-class}
\alias{IdtHetNMANOVA-class}
\alias{IdtLocSNMANOVA-class}
\alias{IdtLocNSNMANOVA-class}
\alias{IdtGenSNMANOVA-class}
\alias{IdtGenNSNMANOVA-class}

\title{Class IdtMANOVA}
\description{IdtMANOVA extends \code{\linkS4class{LRTest}} directly, containing the results of MANOVA tests on the interval-valued data. This class is not used directly, but is the basis for different specializations according to the model assumed for the
distribution in each group. In particular, the following specializations of IdtMANOVA are currently implemented:\cr

\code{\linkS4class{IdtClMANOVA}} extends IdtMANOVA, assuming a classical (i.e., homoscedastic gaussian) setup.\cr

\code{\linkS4class{IdtHetNMANOVA}} extends IdtMANOVA, assuming a heteroscedastic gaussian set-up.\cr

\code{\linkS4class{IdtLocSNMANOVA}} extends IdtMANOVA, assuming a Skew-Normal location model set-up.\cr

\code{\linkS4class{IdtLocNSNMANOVA}} extends IdtMANOVA, assuming either a homoscedastic gaussian or Skew-Normal location model set-up.\cr

\code{\linkS4class{IdtGenSNMANOVA}} extends IdtMANOVA, assuming a Skew-Normal general model set-up.\cr

\code{\linkS4class{IdtGenNSNMANOVA}} extends IdtMANOVA, assuming either a heteroscedastic gaussian or Skew-Normal general model set-up.
}

\section{Slots}{
  \describe{
    \item{\code{NIVar}:}{Number of interval variables.}
    \item{\code{grouping}:}{Factor indicating the group to which each observation belongs to.}
    \item{\code{H0res}:}{Model estimates under the null hypothesis.}
    \item{\code{H1res}:}{Model estimates under the alternative hypothesis.}
    \item{\code{ChiSq}:}{Inherited from class \code{\linkS4class{LRTest}}. Value of the Chi-Square statistics corresponding to the performed test.}
    \item{\code{df}:}{Inherited from class \code{\linkS4class{LRTest}}. Degrees of freedom of the Chi-Square statistics.}
    \item{\code{pvalue}:}{Inherited from class \code{\linkS4class{LRTest}}. p-value of the Chi-Square statistics value, obtained from the Chi-Square distribution with df degrees of freedom.}
    \item{\code{H0logLik}:}{Inherited from class \code{\linkS4class{LRTest}}. Logarithm of the Likelihood function under the null hypothesis.}
    \item{\code{H1logLik}:}{Inherited from class \code{\linkS4class{LRTest}}. Logarithm of the Likelihood function under the alternative hypothesis.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "IdtMANOVA")}: show S4 method for the IdtMANOVA-classes.}
    \item{summary}{\code{signature(object = "IdtMANOVA")}: summary S4 method for the IdtMANOVA-classes.}
    \item{H0res}{\code{signature(object = "IdtMANOVA")}: retrieves the model estimates under the null hypothesis.}
    \item{H1res}{\code{signature(object = "IdtMANOVA")}: retrieves the model estimates under the alternative hypothesis.}
    \item{lda}{\code{signature(x = "IdtClMANOVA")}: Linear Discriminant Analysis using the estimated model parameters.}
    \item{lda}{\code{signature(x = "IdtLocNSNMANOVA")}: Linear Discriminant Analysis using the estimated model parameters.}
    \item{qda}{\code{signature(x = "IdtHetNMANOVA")}: Quadratic Discriminant Analysis using the estimated model parameters.}
    \item{qda}{\code{signature(x = "IdtGenNSNMANOVA")}: Quadratic Discriminant Analysis using the estimated model parameters.}
    \item{snda}{\code{signature(x = "IdtLocNSNMANOVA")}: Discriminant Analysis using maximum likelihood parameter estimates of SkewNormal mixtures assuming a "location" model (i.e., groups differ only in location parameters).}
    \item{snda}{\code{signature(x = "IdtGenSNMANOVA")}: Discriminant Analysis using maximum likelihood parameter estimates of SkewNormal mixtures assuming a general model (i.e., groups differ in all parameters).}
    \item{snda}{\code{signature(x = "IdtGenNSNMANOVA")}: Discriminant Analysis using maximum likelihood parameter estimates of SkewNormal mixtures assuming a general model (i.e., groups differ in all parameters).}
   }
}

\section{Extends}{
Class \code{\linkS4class{LRTest}}, directly.
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Brito, P., Duarte Silva, A. P.  (2012): "Modelling Interval Data with Normal and Skew-Normal Distributions". Journal of Applied Statistics, Volume 39, Issue 1, 3-20.}

\seealso{
  \code{\link{MANOVA}}, \code{\link{lda}}, \code{\link{qda}}, \code{\link{snda}}, \code{\linkS4class{IData}} 
}

\keyword{classes}
\keyword{interval data}
\keyword{MANOVA}
