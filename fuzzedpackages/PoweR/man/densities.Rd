\name{densities} 
\alias{dlaw40}
\alias{dlaw39}
\alias{dlaw38}
\alias{dlaw37}
\alias{dlaw36}
\alias{dlaw35}
\alias{dlaw34}
\alias{dlaw33}
\alias{dlaw32}
\alias{dlaw31}
\alias{dlaw30}
\alias{dlaw29}
\alias{dlaw28}
\alias{dlaw27}
\alias{dlaw26}
\alias{dlaw25}
\alias{dlaw24}
\alias{dlaw23}
\alias{dlaw22}
\alias{dlaw21}
\alias{dlaw20}
\alias{dlaw19}
\alias{dlaw18}
\alias{dlaw17}
\alias{dlaw16}
\alias{dlaw15}
\alias{dlaw14}
\alias{dlaw13}
\alias{dlaw12}
\alias{dlaw11}
\alias{dlaw10}
\alias{dlaw9}
\alias{dlaw8}
\alias{dlaw7}
\alias{dlaw6}
\alias{dlaw5}
\alias{dlaw4}
\alias{dlaw3}
\alias{dlaw2}
\alias{dlaw1}
\alias{indicator}

\title{Density function}

\description{Evaluate the density function at a vector points.}

\details{
  Use the function by typing:
  
  dlawj(x,par1,par2,etc.)

  where \eqn{j} is the index of the law and par1, par2, etc. are the
  parameters of law \eqn{j}.

  The \code{indicator} function takes a vector \code{x} of length \eqn{n} as first
  argument and two real values \eqn{a<b}. It returns a vector of length \eqn{n}
  which contains only 0s and 1s (1 if the corresponding value in \code{x} is
  strictly between \code{a} and \code{b}).
}


\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical
Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03
}
\author{P. Lafaye de Micheaux, V. A. Tran
}
\keyword{distributions}
