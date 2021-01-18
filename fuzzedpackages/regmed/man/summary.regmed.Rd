 \name{summary.regmed}
\alias{summary.regmed}
\title{
Summary of regmed object
}
\description{
Summarizes alpha and beta coefficients of mediation effect, and delta
for direct effect, in a regularized mediation model (regmed object). 
}
\usage{
\method{summary}{regmed}(object, \dots)
}
\arguments{
  \item{object}{ regmed object}
  \item{\dots}{ optional arguments for summary method}
}
\details{
Prints summary of model fit. Alpha and beta represent the indirect mediation
effects:  x[alpha] -> mediator -> [beta]outcome. Delta represents the
direct effect: x[delta] -> y. Var(x) and Var(y) represent the residual
variances. The sum of alpha*beta represents the total mediation effect.
}
\value{ invisible, data.frame for alpha/beta pairs in fitted object
}
\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}
