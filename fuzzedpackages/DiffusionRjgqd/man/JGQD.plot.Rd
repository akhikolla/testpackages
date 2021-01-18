\name{JGQD.plot}
\alias{JGQD.plot}
\title{
Quick Plots for DiffusionRjgqd Objects
}
\description{
\code{JGQD.plot()} recognizes output objects calculated using routines from the \bold{DiffusionRjgqd} package and subsequently constructs an appropriate plot, for example a perspective plot of a transition density.
}
\usage{
JGQD.plot(x, thin = 1, burns, h = FALSE,  palette = 'mono')
}
\arguments{
  \item{x}{
Generic JGQD-objects, i.e. \code{res = JGQD.density()}.
}
  \item{thin}{
Thinning interval for \code{.mcmc} objects.
}
  \item{burns}{
Number of parameter draws to discard for \code{.mcmc} objects.
}
  \item{h}{
  if \code{TRUE} a histogram is drawn i.s.o. a trace plot.
}
\item{palette}{Colour palette for drawing trace plots. Default \code{palette = 'mono'}, otherwise a qualitative palette will be used.}
}

\value{Varies in accordance with input type.
}
\author{
Etienne A.D. Pienaar: \email{etiannead@gmail.com}
}
\references{
Updates available on GitHub at \url{https://github.com/eta21}.
}


\seealso{
\code{\link{JGQD.mcmc}}, \code{\link{JGQD.density}}, \code{\link{BiJGQD.density}} etc.
}
\examples{
  \donttest{
#===============================================================================
# Plot the transitional density of a jump diffusion
#-------------------------------------------------------------------------------
rm(list=ls(all=TRUE))
library(DiffusionRjgqd)

JGQD.remove()
# Define the jump diffusion using the DiffusionRjgqd syntax:
G1=function(t){0.2*5+0.1*sin(2*pi*t)}
G2=function(t){-0.2}
Q1=function(t){0.2}

# State dependent intensity:
Lam0 = function(t){1}
Lam1    = function(t){0.1}

# Normally distributed jumps: N(1,0.2)
Jmu    = function(t){1.0}
Jsig   = function(t){0.2}
# Normal distribution is the default:
res_1  = JGQD.density(4,seq(2,10,1/10),0,2.5,1/100,factorize=FALSE)

JGQD.plot(res_1)

  }
}
\keyword{plot}

