\name{JGQD.mcmc}
\alias{JGQD.mcmc}

\title{
MCMC Inference on Jump Generalized Quadratic Diffusion Models (JGQDs).
}
\description{
\code{JGQD.mcmc()} uses parametrised coefficients (provided by the user as R-functions) to construct a C++ program in real time that allows the user to perform Bayesian inference on the resulting jump diffusion model. Given a set of starting parameters, a MCMC chain is returned for further analysis.
The structure of the model is predefined and coefficients may be provided for models nested within the generalized quadratic diffusion framework.
\if{html}{
 \code{JGQD.mcmc()} performs inference using the Metropolis-Hastings algorithm for jump diffusions of the form:

 \figure{ScalarEqn1.png}

 where

 \figure{ScalarEqn2.png}

  \figure{ScalarEqn3.png}

  and

   \figure{ScalarEqn4.png}

  describes a Poisson process with jumps of the form:

   \figure{ScalarEqn6.png}

  arriving with intensity

   \figure{ScalarEqn5.png}

  subject to a jump distribition of the form:

  \figure{ScalarEqn7.png}

 }
\if{latex}{
\code{JGQD.mcmc()} performs inference using the Metropolis-Hastings algorithm for jump diffusions of the form:
\deqn{dX_t = \mu(X,t)dt +\sigma(X_t,t)dW_t+dP_t\\}
where
\deqn{\mu(X,t) = G_0(t,\theta)+G_1(t,\theta)X_t+G_2(t,\theta)X_t^2,\\}
\deqn{\sigma(X,t)=\sqrt{Q_0(t,\theta)+Q_1(t,\theta)X_t+Q_2(t,\theta)X_t^2,}\\}
and \eqn{dP_t=j(X_t,\dot{z}_t)dN_t} describes a Poisson process with intensity:
\deqn{\lambda(X,t)={\lambda_0(t,\theta)+\lambda_1(t,\theta)X_t+\lambda_2(t,\theta)X_t^2},\\}
jump coefficient
\deqn{j(X_t,\dot{z}_t) = \dot{z}_t} or \deqn{j(X_t,\dot{z}_t) = \dot{z}_tX_t\\ }
and
\deqn{\dot{z}_t\sim\mbox{Normal}(\mu_J(t,\theta),\sigma_J^2(t,\theta)),\\}
\deqn{\dot{z}_t\sim\mbox{Exponential}(\lambda_J(t,\theta)),\\}
\deqn{\dot{z}_t\sim\mbox{Gamma}(\alpha_J(t,\theta),\beta_J(t,\theta)),\\}
or
\deqn{\dot{z}_t\sim\mbox{Laplace}(a_J(t,\theta),b_J(t,\theta)).\\}
 }

}
\usage{
JGQD.mcmc(X, time, mesh = 10, theta, sds, updates = 1000,
          burns = min(round(updates/2), 25000), Jtype = "Add", Jdist = "Normal",
          Dtype = "Saddlepoint", RK.order = 4, exclude = NULL, plot.chain = TRUE,
          wrt = FALSE, Tag = NA, factorize = TRUE, print.output = TRUE,
          decode = TRUE, palette = 'mono')
}

\arguments{
  \item{X}{
Time series (vector) of discretely observed points of the process of interest. These may be non-equidistant observations (see \code{time}).
}
  \item{time}{
A vector of time-stamps associated with each observation in \code{X}.
}
  \item{mesh}{
The number mesh points between any two given data points.
}
  \item{theta}{
The parameter vector of the process. \code{theta} are taken as the starting values of the MCMC chain and gives the dimension of the parameter vector used to calculate the DIC. Care should be taken to ensure that each element in \code{theta} is in fact used within the coefficient-functions, otherwise redundant parameters will be counted in the calculation of the DIC.
}
  \item{sds}{
  Proposal distribution standard deviations. That is, for the i-th parameter the proposal distribution is \emph{~ Normal(...,\code{sds}[i]^2)}
}
  \item{updates}{
The number of chain updates (including burned updates) to perform.
}
  \item{burns}{
The number of updates to burn. That is, the first \code{burns} values are omitted from the inference, although the entire chain is returned.
}
  \item{exclude}{ Vector indicating which transitions to exclude from the analysis. Default = \code{NULL}.}
  \item{plot.chain}{
  If \code{TRUE} (default), a trace plot is made of the resulting MCMC chain (see details).
}
\item{RK.order}{The order of the Runge-Kutta solver used to approximate the trajectories of cumulants. Must be 4 or (default) 10.}
  \item{Dtype}{Character string indicating the type of density approximation (see details) to use. Types: \code{'Saddlepoint'} is supported in the current version of the software.}
\item{Tag}{\code{Tag} can be used to name (tag) an MCMC run e.g. \code{Tag='Run_1'}}
\item{wrt}{If \code{TRUE} a .cpp file will be written to the current directory. For bug report diagnostics.}
\item{Jdist}{Valid entries are 'Normal', 'Expnential', 'Gamma' and 'Laplace'.}
\item{Jtype}{Valid types are 'Add' or 'Mult'.}
\item{factorize}{Should factorization be used (default = TRUE).}
\item{print.output}{If \code{TRUE}, model information is printed to the console.}
\item{decode}{Should the algorithm estimate jump detection probabilities? Default value is \code{TRUE}. }
\item{palette}{Colour palette for drawing trace plots. Default \code{palette = 'mono'}, otherwise a qualitative palette will be used.}
}
\details{

\code{JGQD.mcmc()} operates by searching the workspace for functions with names that match the coefficients of the predefined stochastic differential equation. Only the required coefficients need to be specified e.g. \strong{\code{G0(t)}},\strong{\code{G1(t)}} and \strong{\code{Q0(t)}} for an Ornstein-Uhlenbeck model. Unspecified coefficients are ignored. When a new model is to be defined, the current model may be removed from the workspace by using the \code{\link{JGQD.remove}} function, after which the new coefficients may be supplied.
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
Etienne A.D. Pienaar: \email{etiennead@gmail.com}
}
\note{
\bold{Note [1]:} When \code{plot.chain} is \code{TRUE}, a trace plot is created of the resulting MCMC along with the acceptance rate at each update. This may save time when scrutinizing initial MCMC runs.
}

\section{Syntactical jargon}{
\bold{Synt. [1]:} The coefficients of the JGQD may be parameterized using the reserved variable \code{theta}. For example:

\code{G0 <- function(t){theta[1]*(theta[2]+sin(2*pi*(t-theta[3])))}}.


\bold{Synt. [2]:} Due to syntactical differences between R and C++ special functions have to be used when terms that depend on \code{t}. When the function cannot be separated in to terms that contain a single \code{t}, the \code{prod(a,b)} function must be used. For example:

\code{G0 <- function(t){0.1*(10+0.2*sin(2*pi*t)+0.3*prod(sqrt(t),1+cos(3*pi*t)))}}.

Here sqrt(t)*cos(3*pi*t) constitutes the product of two terms that cannot be written i.t.o. a single \code{t}. To circumvent this isue, one may use the \code{prod(a,b)} function.

\bold{Synt. [3]:} Similarly, the ^ - operator is not overloaded in C++. Instead the \code{pow(x,p)} function may be used to calculate x^p. For example sin(2*pi*t)^3 in:

\code{G0 <- function(t){0.1*(10+0.2*pow(sin(2*pi*t),3))}}.

}

\seealso{
\code{\link{JGQD.remove}},  \code{\link{BiJGQD.mcmc}}.
}


\examples{
\donttest{

#===============================================================================

 # For detailed notes and examples on how to use the JGQD.mcmc() function, see
 # the following vignette:

 RShowDoc('Part_4_Likelihood_Inference',type='html','DiffusionRjgqd')

#===============================================================================
}
}
\keyword{syntax}
\keyword{C++}
\keyword{mcmc}


