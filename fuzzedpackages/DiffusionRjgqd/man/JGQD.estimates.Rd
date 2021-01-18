\name{JGQD.estimates}
\alias{JGQD.estimates}
\title{
Extract Parmaeter Estimates from \code{.mle()} or \code{.mcmc()} Objects.
}
\description{
\code{JGQD.estimates()} calculates parameter estimates from \code{.mle()} or \code{.mcmc()} model objects.
}
\usage{
JGQD.estimates(x, thin = 100, burns, CI = c(0.05, 0.95), corrmat =
                 FALSE, acf.plot = TRUE, palette='mono')
}

\arguments{
  \item{x}{List object of type 'JGQD.mcmc' or 'JGQD.mle'. That is, when \code{model =JGQD.mcmc()} then \code{model} constitutes an appropriate object for \code{x}.
}
  \item{thin}{
Thinnging level for parameter chain.
}
  \item{burns}{
Number of MCMC updates to discard before calculating estimates.
}
 \item{CI}{
Credibility interval quantiles (for MCMC chains).
}

  \item{corrmat}{
If TRUE, an estimated correlation matrix is returned in addition to the estimate vector.
}
  \item{acf.plot}{
If TRUE, an acf plot is drawn for each element of the parameter chain.
}
 \item{palette}{Colour palette for drawing trace plots. Default \code{palette = 'mono'}, otherwise a qualitative palette will be used.}
}

\value{
Data frame with parameter estimates and appropriate interval statistics.
}

\author{
Etienne A.D. Pienaar: \email{etiannead@gmail.com}
}
\references{
Updates available on GitHub at \url{https://github.com/eta21}.
}
\examples{
\donttest{
example(JGQD.mcmc)
}
}


\seealso{
\code{\link{JGQD.mcmc}},  \code{\link{BiJGQD.mcmc}}.
}

