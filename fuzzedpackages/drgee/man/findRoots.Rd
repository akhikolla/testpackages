\encoding{latin1}
\name{findRoots}
\alias{findRoots}

\title{
Non Linear Equation System Solving
}
\description{
	A wrapper around \code{nleqslv} from the \pkg{nleqslv} package
	to solve a non linear system of equations.
}
\usage{
findRoots(beta.init, eq.func, d.eq.func = NULL, arg.list, ...)
}
\arguments{
  \item{beta.init}{
    An initial guess for the zero.
}
  \item{eq.func}{
    A function of two variables for which the zero are sought. Its first
    argument \code{beta} should be a vector over which the zeros are sought and
    the second argument \code{arg.list} a list of additional arguments.
}
  \item{d.eq.func}{
    A function to return the Jacobian of \code{eq.func} taking the same
    arguments as \code{eq.func}. Supplying this function can speed up
    calculations. Default is \code{NULL}.
}
  \item{arg.list}{
	The second argument to \code{eq.func} and \code{d.eq.func}}
  \item{...}{
	A list of additional arguments to be passed to \code{nleqslv}
	}
}
\details{
  \code{findRoots} calculates zeros fo the function \code{eq.func}
  and is the default equation solving function in \code{drgee}.
  It is supplied as a separate function in order to allow users
  to use other equation solvers by writing their own wrapper
  with the same interface as \code{findRoots}.
}
\value{
	The value is a list containing the following arguments:
	\item{roots }{The zero(s) of the function \code{eq.func}.
	}
	\item{optim.object }{The optimization object returned from \code{nleqslv}.
	}
}
\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}
\seealso{
\code{\link{nleqslv}} in package \pkg{nleqslv}.
}
\keyword{optimize}
