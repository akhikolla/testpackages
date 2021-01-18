\name{JSDEsim1}
\alias{JSDEsim1}
\docType{data}
\title{
Simulated Trajectory of a Scalar Jump Diffusion.
}
\description{
Simulated trajectory of a scalar jump diffusion with Normal distributed jumps ariving with stochastic intensity.
\if{latex}{The dynamics of the process is given by the jump SDE:
\deqn{dX_t = (5-X_t)dt+0.1\sqrt{X_t}dB_t+dP_t\\}
where \eqn{dP_t=\dot{z}_tdN_t} describes a Poisson process with intensity:
\deqn{\lambda(X_t) = 0.5X_t\\}
and
\deqn{\dot{z}\sim\mbox{N}(0.5,0.5^2)\\}}

\if{html}{
 The dynamics of the process is given by the jump SDE:

 \figure{JSDEsim1.png}{jeqn}
}

}
\usage{data("JSDEsim1")}
\format{
  A data frame with 251 observations on the following 2 variables.
  \describe{
    \item{\code{Xt}}{Observed trajectory.}
    \item{\code{time}}{Time index for the observed trajectory}
    \item{\code{contains_jump}}{Indicator variable to indicate if a given transition contains a jump (ignore the first value).}
  }
}

\examples{
data(JSDEsim1)
}
\keyword{datasets}
