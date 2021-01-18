\name{boost.mutex-class}
\Rdversion{1.1}
\docType{class}
\alias{mutex-class}
\alias{boost.mutex-class}
\alias{is.timed,boost.mutex-method}
\alias{lock.shared,boost.mutex-method}
\alias{lock,boost.mutex-method}
\alias{shared.name,boost.mutex-method}
\alias{timeout,boost.mutex-method}
\alias{unlock,boost.mutex-method}

\title{Class ``boost.mutex''}
\description{The \code{boost.mutex} class provides an \R interface
to the mutex functionality implemented in the Boost C++ library}
\section{Objects from the Class}{
Unlike many \R objects, objects should not be created by calls of the form
\code{new("boost.mutex", ...)}.  The function \code{boost.mutex()}
is intended for the user.
}
\section{Slots}{
	 \describe{
    \item{\code{isRead}:}{This is used internally to maintain state information
                          and should not be touched by a user.}
    \item{\code{mutexInfoAddr}:}{Object of class \code{"externalptr"}
which keeps track of information relevant to the mutex.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{mutex}"}, directly.
}
\section{Methods}{
  \describe{
    \item{describe}{\code{signature(x = "boost.mutex")}: ... }
    \item{is.timed}{\code{signature(m = "boost.mutex")}: ... }
    \item{lock.shared}{\code{signature(m = "boost.mutex")}: ... }
    \item{lock}{\code{signature(m = "boost.mutex")}: ... }
    \item{shared.name}{\code{signature(m = "boost.mutex")}: ... }
    \item{timeout}{\code{signature(m = "boost.mutex")}: ... }
    \item{unlock}{\code{signature(m = "boost.mutex")}: ... }
	 }
}
\author{Michael J. Kane <bigmemoryauthors@gmail.com>}
\seealso{\code{\link{boost.mutex}}}
\examples{
showClass("boost.mutex")
}
\keyword{classes}
