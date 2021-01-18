\name{Classes of basic handling of genealogy data}
\Rdversion{1.1}
\docType{class}
\alias{GLmultiVector-class}
\alias{GLmultiArray4-class}
\alias{GLmultiMatrix-class}
\alias{GLmultiNumber-class}
\alias{[<-,GLmultiVector,ANY,ANY,ANY-method}
\alias{[,GLmultiVector,ANY,ANY,ANY-method}
\alias{depth,GLmultiVector-method}
\alias{Dim,GLmultiVector-method}
\alias{[<-,GLmultiArray4,ANY,ANY,ANY-method}
\alias{[,GLmultiArray4,ANY,ANY,ANY-method}
\alias{depth,GLmultiArray4-method}
\alias{Dim,GLmultiArray4-method}
\alias{[<-,GLmultiMatrix,ANY,ANY,ANY-method}
\alias{[,GLmultiMatrix,ANY,ANY,ANY-method}
\alias{depth,GLmultiMatrix-method}
\alias{Dim,GLmultiMatrix-method}
\alias{[<-,GLmultiNumber,ANY,ANY,ANY-method}
\alias{[,GLmultiNumber,ANY,ANY,ANY-method}
\alias{depth,GLmultiNumber-method}
\alias{Dim,GLmultiNumber-method}

\title{Classes \code{"GLmultiVector"} \code{"GLmultiArray4"} \code{"GLmultiMatrix"} \code{"GLmultiNumber"}}
\description{
Objects created to carry information mostly between invisible functions.
}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("GLmultiVector", 'matrix', depth = 'integer')}\cr
\code{new("GLmultiArray4", 'Array', depth = 'integer')}\cr
\code{new("GLmultiMatrix", 'Array', depth = 'integer')}\cr
\code{new("GLmultiNumber", 'numeric', depth = 'integer', .Names = 'character')}
}
\section{Slots}{
  \describe{
    GLmultiVector
    \item{\code{.Data}:}{Object of class \code{"matrix"} }
    \item{\code{depth}:}{Object of class \code{"integer"} }\cr
    GLmultiArray4, GLmultiMatrix
    \item{\code{.Data}:}{Object of class \code{"array"} }
    \item{\code{depth}:}{Object of class \code{"integer"} }\cr
    GLmultiNumber
    \item{\code{.Data}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{depth}:}{Object of class \code{"integer"} ~~ }
    \item{\code{.Names}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "GLmultiVector", i = "ANY", j = "ANY", value = "ANY")}: ... }
    \item{[}{\code{signature(x = "GLmultiVector", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{depth}{\code{signature(x = "GLmultiVector")}: ... }
    \item{Dim}{\code{signature(object = "GLmultiVector")}: ... }
	 }
}

\examples{
showClass("GLmultiVector")
showClass("GLmultiArray4")
showClass("GLmultiMatrix")
showClass("GLmultiNumber")
}

\keyword{classes}

