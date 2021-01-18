\name{JGQD.remove}
\alias{JGQD.remove}

\title{
Remove the Coefficients of a JGQD Model.
}
\description{Removes any existing coefficient functions from the current workspace.
}
\usage{
JGQD.remove()
}

\details{
\code{JGQD.remove} clears the workspace of functions with names that match the coefficients of the 1D JGQD. This may be used when more than one model is specified  in a given session.
}
\value{
No value is returned.
}
\references{
Updates available on GitHub at \url{https://github.com/eta21}.
}
\author{
Etienne A.D. Pienaar: \email{etiennead@gmail.com}
}
\note{
\code{JGQD.remove} simply searches the workspace for functions with definitions that match the form of the DiffusionRjgqd interface and removes them from the workspace, freeing up the user to redefine a diffusion with new coefficients.
}



\seealso{
\code{\link{JGQD.density}} and \code{\link{BiJGQD.density}}.
}
\keyword{remove models}


