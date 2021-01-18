\name{GLgroup-class}
\docType{class}
\alias{GLgroup-class}
\alias{initialize,GLgroup-method}
\alias{[,GLgroup,ANY,ANY,ANY-method}

\title{Class \code{"GLgroup"}}
\description{
 Object representing a set of proband in different groups.\cr
 Each element of this list is a group named with the proband number.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GLgroup", 'list')}.
}
\section{Extends}{
Class \code{"\linkS4class{GLgen}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Data=list)}: ... }
    \item{[}{\code{signature(object = "GLgroup",ANY,ANY,ANY)}:... }
 }
}

\examples{
showClass("GLgroup")
}
\keyword{classes}
