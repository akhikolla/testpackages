\name{BiJGQD.mcmc}
\alias{BiJGQD.mcmc}

\title{
MCMC Inference on Bivariate Jump Generalized Quadratic Diffusions (2D JGQDs).
}
\description{
\code{BiJGQD.mcmc()} uses parametrised coefficients (provided by the user as R-functions) to construct a C++ program in real time that allows the user to perform Bayesian inference on the resulting diffusion model. Given a set of starting parameters and other input parameters, a MCMC chain is returned for further analysis. The user may specify any model within the JGQD framework by defining parametrised functions giving the form of the coefficients of the model.
 \if{html}{
That is, \code{BiJGQD.density} generates approximate transitional densities for a class of bivariate jump diffusion processes with SDE:

\figure{BivJumpDiff1.png}

where

\figure{BivJumpDiff4.png}

\figure{BivJumpDiff3.png}

and

\figure{BivJumpDiff2.png}

describes a bivariate Poisson process with jump matrix:

\figure{BivJumpDiff5.png}

and intensity vector

\figure{BivJumpDiff6.png}

}
\if{latex}{
\code{BiJGQD.density()} approximates the transitional density of jump diffusions of the form:
\deqn{d\mathbf{X}_t = \boldsymbol\mu(\mathbf{X}_t,t)dt +\boldsymbol\sigma(\mathbf{X}_t,t)d\mathbf{W}_t+d\mathbf{P}_t\\}
where
\deqn{\mu^{(1)}(\mathbf{X}_t,t)= \sum_{i+j\leq 2}a_{ij}(t)X_t^iY_t^j,\\}
\deqn{\mu^{(2)}(\mathbf{X}_t,t)= \sum_{i+j\leq 2}b_{ij}(t)X_t^iY_t^j,\\}
and
\deqn{\boldsymbol\sigma(X_t,Y_t,t)\boldsymbol\sigma^{\prime}(X_t,Y_t,t) = (\gamma^{(i,j)}(\mathbf{X}_t,t))_{i,j=1,2}\\}
such that
\deqn{\gamma^{(1,1)}(\mathbf{X}_t,t)= \sum_{i+j\leq 2}c_{ij}(t)X_t^iY_t^j,\\}
\deqn{\gamma^{(1,2)}(\mathbf{X}_t,t)= \sum_{i+j\leq 2}d_{ij}(t)X_t^iY_t^j,\\}
\deqn{\gamma^{(2,1)}(\mathbf{X}_t,t)= \sum_{i+j\leq 2}e_{ij}(t)X_t^iY_t^j,\\}
\deqn{\gamma^{(2,2)}(\mathbf{X}_t,t)= \sum_{i+j\leq 2}f_{ij}(t)X_t^iY_t^j,\\}
and \eqn{d\mathbf{P}_t=J(\mathbf{X}_t,t)d\mathbf{N}_t} describes a Poisson process with jump matrix:
\deqn{\mathbf{J}(\mathbf{X}_t,t,\dot{\mathbf{z}}_t)^{(1,1)}=\dot{z}_{11},\\}
\deqn{\mathbf{J}(\mathbf{X}_t,t,\dot{\mathbf{z}}_t)^{(2,1)}=\dot{z}_{21},\\}
if jumps are additive or
\deqn{\mathbf{J}(\mathbf{X}_t,t,\dot{\mathbf{z}}_t)^{(1,1)}=\dot{z}_{11}X_t,\\}
\deqn{\mathbf{J}(\mathbf{X}_t,t,\dot{\mathbf{z}}_t)^{(2,1)}=\dot{z}_{21}Y_t,\\}
if jumps are multiplicative. For purposes of this package the arrival rate is decomposed as:
\deqn{\lambda(\mathbf{X}_t,t)= \sum_{i+j\leq 1}\lambda_{ij}(t)X_t^iY_t^j\\}
and
 \deqn{(\dot{z}_{11},\dot{z}_{21})^\prime \sim \mbox{MVN}(\boldsymbol\mu_J,\boldsymbol\Sigma_J).\\}
 }
}
\usage{
BiJGQD.mcmc(X, time, mesh = 10, theta, sds, updates = 10,
            burns = min(round(updates/2),25000), exclude = NULL, plot.chain = TRUE,
            RK.order = 4, wrt = FALSE, Tag = NA, Dtype = "Saddlepoint",
            Jdist = "MVNormal", Jtype ='Add', adapt = 0, print.output = TRUE,
            decode = TRUE, palette = 'mono')
}
\arguments{
  \item{X}{A matrix containing rows of data points to be modelled. Although observations are allowed to be non-equidistant, observations in both dimensions are assumed to occur at the same time epochs (i.e. \code{time} gives the time signature for both dimensions).
}
  \item{time}{A vector containing the time epochs at which observations were made.
}
  \item{mesh}{The number of mesh points in the time discretization.
}
  \item{theta}{The parameter vector of the process. theta are taken as the starting values of the MCMC chain and gives the dimension of the parameter vector used to calculate the DIC. Care should be taken to ensure that each element in theta is in fact used within the coefficient-functions, otherwise redundant parameters will be counted in the calculation of the DIC.
}
  \item{sds}{ Proposal distribution standard deviations. That is, for the i-th parameter the proposal distribution is \emph{~ Normal(...,\code{sds}[i]^2).
}}
  \item{updates}{The number of MCMC updates/iterations to perform (including burn-in).
}
  \item{burns}{The number of MCMC updates/iterations to burn.}
  \item{exclude}{ Vector indicating which transitions to exclude from the analysis. Default = \code{NULL}.}
  \item{plot.chain}{If = TRUE (default), a trace plot of the MCMC chain will be made along with a trace of the acceptance rate.}
\item{RK.order}{The order of the Runge-Kutta solver used to approximate the trajectories of cumulants. Must be 4 (default) or 10.}
\item{Tag}{\code{Tag} can be used to name (tag) an MCMC run e.g. \code{Tag='Run_1'}}
\item{Dtype}{The density approximant to use. Valid types are \code{"Saddlepoint"} (default), \code{"Edgeworth"} or \code{"Normal"}.}
\item{Jdist}{Valid entries are 'MVNormal' (currently).}
\item{Jtype}{Valid types are 1 or 2.}
\item{adapt}{For development purposes.}
\item{wrt}{If \code{TRUE} a .cpp file will be written to the current directory. For bug report diagnostics.}
\item{print.output}{If \code{TRUE} information about the model and algorithm is printed to the console. }
\item{decode}{Should the algorithm estimate jump detection probabilities? Default value is \code{TRUE}. }
\item{palette}{Colour palette for drawing trace plots. Default \code{palette = 'mono'}, otherwise a qualitative palette will be used.}
}
\value{
  \item{par.matrix }{A matrix containing the MCMC chain on \code{theta}. }
  \item{acceptence.rate }{A vector containing the acceptance rate of the MCMC at every iteration. }
  \item{model.info }{A list of variables pertaining to inference calculations.}
  \item{model.info$elapsed.time}{The runtime, in h/m/s format,of the MCMC procedure (excluding compile time).}
  \item{model.info$time.homogeneous}{`No' if the model has time-homogeneous coefficients and `Yes' otherwise.}
  \item{model.info$p}{The dimension of \code{theta}.}
  \item{model.info$DIC}{Calculated Deviance Information Criterion.}
  \item{model.info$pd}{Effective number of parameters (see \code{model.info$DIC}).}
  \item{decode.prob}{Estimated jump detection probabilities.}
}
\references{
Updates available on GitHub at \url{https://github.com/eta21}.

Daniels, H.E. 1954 Saddlepoint approximations in statistics. \emph{Ann. Math. Stat.}, \bold{25}:631--650.

Eddelbuettel, D. and Romain, F. 2011 Rcpp: Seamless R and C++ integration. \emph{Journal of Statistical Software}, \bold{40}(8):1--18,. URL \url{http://www.jstatsoft.org/v40/i08/}.

Eddelbuettel, D. 2013 \emph{Seamless R and C++ Integration with Rcpp}. New York: Springer. ISBN
978-1-4614-6867-7.

Eddelbuettel, D. and Sanderson, C. 2014 Rcpparmadillo: Accelerating R with high-performance C++
linear algebra. \emph{Computational Statistics and Data Analysis}, \bold{71}:1054--1063. URL
\url{http://dx.doi.org/10.1016/j.csda.2013.02.005}.

Feagin, T. 2007 A tenth-order Runge-Kutta method with error estimate. \emph{In Proceedings of the IAENG
Conf. on Scientifc Computing}.

Varughese, M.M. 2013 Parameter estimation for multivariate diffusion systems. \emph{Comput. Stat. Data An.},
\bold{57}:417--428.
}
\author{
Etienne A.D. Pienaar \email{etiannead@gmail.com}
}

\note{

\bold{Note [1]:} When \code{plot.chain} is \code{TRUE}, a trace plot is created of the resulting MCMC along with the acceptance rate at each update. This may save time when
scrutinizing initial MCMC runs.

}

\section{Syntactical jargon}{
\bold{Synt. [1]:} The coefficients of the 2D JGQD may be parameterized using the reserved variable \code{theta}. For example:

\code{a00 <- function(t){theta[1]*(theta[2]+sin(2*pi*(t-theta[3])))}}.


\bold{Synt. [2]:} Due to syntactical differences between R and C++ special functions have to be used when terms that depend on \code{t}. When the function cannot be separated in to terms that contain a single \code{t}, the \code{prod(a,b)} function must be used. For example:

\code{a00 <- function(t){0.1*(10+0.2*sin(2*pi*t)+0.3*prod(sqrt(t),1+cos(3*pi*t)))}}.

Here sqrt(t)*cos(3*pi*t) constitutes the product of two terms that cannot be written i.t.o. a single \code{t}. To circumvent this isue, one may use the \code{prod(a,b)} function.

\bold{Synt. [3]:} Similarly, the ^ - operator is not overloaded in C++. Instead the \code{pow(x,p)} function may be used to calculate x^p. For example sin(2*pi*t)^3 in:

\code{a00 <- function(t){0.1*(10+0.2*pow(sin(2*pi*t),3))}}.
}

\section{Warning }{

\bold{Warning [1]}: The parameter \code{mesh} is used to discretize the transition horizons between successive observations. It is thus important to ensure that \code{mesh} is not too small when large time differences are present in \code{time}. Check output for max(dt) and divide by \code{mesh}.

}

\seealso{

\code{\link{JGQD.remove}} and \code{\link{JGQD.mcmc}}.

}

\examples{
\donttest{

#===============================================================================

 # For detailed notes and examples on how to use the BiJGQD.mcmc() function, see
 # the following vignette:

 RShowDoc('Part_4_Likelihood_Inference',type='html','DiffusionRjgqd')

#===============================================================================
 }
 }
\keyword{Syntax}
\keyword{C++}
\keyword{MCMC}
