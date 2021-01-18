\name{GLgen-class}
\docType{class}
\alias{GLgen-class}
\alias{initialize,GLgen-method}
\alias{depth,GLgen-method}
\alias{length,GLgen-method}

\title{Class \code{"GLgen"}}
\description{
 Object containing a genealogy.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{gen.genealogy('data.frame')}.
}
\section{Extends}{
Class \code{"\linkS4class{GLgroup}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Data="integer",Date="character")}: ... }
    \item{depth}{\code{signature(x = "GLgen")}}{ Gives the depth of the genealogy. }
    \item{length}{\code{signature(x = "GLgen")}}{ Gives the number of individuals in the genealogy. }
  }
}

\examples{
showClass("GLgen")
}
\keyword{classes}
