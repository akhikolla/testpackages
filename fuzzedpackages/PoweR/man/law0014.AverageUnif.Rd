\name{law0014.AverageUnif}
\alias{law0014.AverageUnif}
\alias{law0014}

\title{The Average Uniform Distribution}

\description{Random generation for the Average Uniform distribution with parameters \code{size}, \code{a} and \code{b}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{size}, \code{a} and \code{b} are not specified they assume the default values of 2, 0 and 1.

The Average Uniform distribution has density:
\deqn{ \frac{k^k}{(k-1)!}\sum_{j=0}^{\lfloor k\frac{x-a}{b-a} \rfloor}(-1)^j{k \choose j}(\frac{x-a}{b-a}-\frac{j}{k})^{k-1} }{ 
       f(x) = k^k/(k-1)! \sum_{j=0}^{\lfloor k (x-a)/(b-a) \rfloor} (-1)^j choose(k,j) ((x-a)/(b-a) - j/k)^(k-1) }
where \code{size = k} and for \eqn{a \le x \le b}.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 

Quesenberry and Miller (1977), Power studies of some tests for uniformity, \emph{Journal of Statistical Computation and Simulation}, \bold{5}(3), 169--191 (see p. 179)
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{\code{\link{law0007.Uniform}} for the Uniform distribution.

\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(14,10000,law.pars=c(9,2,3))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Average Uniform}
