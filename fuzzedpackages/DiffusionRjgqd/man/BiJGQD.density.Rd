\name{BiJGQD.density}
\alias{BiJGQD.density}
\title{
Generate the Transition Density of a Bivariate Jump Generalized Quadratic Diffusion Model (2D JGQD).
}
\description{

\code{BiJGQD.density} generates approximate transitional densities for bivariate generalized quadratic jump diffusions (JGQDs). Given a starting coordinate, (\code{Xs}, \code{Ys}), the approximation is evaluated over a lattice \code{Xt} x \code{Yt} for an equispaced discretization (intervals of width \code{delt}) of the transition time horizon [\code{s}, \code{t}].
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
BiJGQD.density(Xs, Ys, Xt, Yt, s, t, delt, Dtype, Jdist, Jtype, print.output,
               eval.density)
}
\arguments{
  \item{Xt}{x-Coordinates of the lattice at which to evaluate the transition density.}
  \item{Yt}{y-Coordinates of the lattice at which to evaluate the transition density.}
  \item{Xs}{Initial x-coordinate.}
  \item{Ys}{Initial y-coordinate.}
  \item{s}{Starting time of the diffusion.}
  \item{t}{Final time at which to evaluate the transition density.}
  \item{delt}{Step size for numerical solution of the cumulant system. Also used for the discretization of the transition time-horizon. See warnings \bold{[1]} and \bold{[2]}.}
  \item{Dtype}{The density approximant to use. Valid types are \code{"Saddlepoint"} (default) or \code{"Edgeworth"}.}
  \item{Jdist}{Valid entries are 'MVNormal' (currently).}
  \item{Jtype}{Valid types are 1 or 2.}
  \item{print.output}{If \code{TRUE} information about the model and algorithm is printed to the console. }
 \item{eval.density}{If \code{TRUE}, the density is evaluated in addition to calculating the moment eqns.}
}
\details{


\if{html}{\figure{TransDens2.png}{GQD}}

\code{BiJGQD.density} constructs an approximate transition density for a class of quadratic diffusion models. This is done by first evaluating the trajectory of the cumulants/moments of the diffusion numerically as the solution of a system of ordinary differential equations over a time horizon \code{[s,t]} split into equi-distant points \code{delt} units apart. Subsequently, the resulting cumulants/moments are carried into a density approximant (by default, a saddlepoint approximation) in order to evaluate the transtion surface.
}
\value{
  \item{density}{3D Array containing approximate density values. Note that the 3rd dimension represents time.}
  \item{Xt}{Copy of x-coordinates.}
  \item{Yt}{Copy of y-coordinates.}
  \item{time}{A vector containing the time mesh at which the density was evaluated.}
  \item{cumulants}{A matrix giving the cumulants of the diffusion. Cumulants are indicated by row-names.}
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
Etienne A.D. Pienaar: \email{etiannead@gmail.com}
}


\section{Warning}{
\bold{Warning [1]}:
The system of ODEs that dictate the evolution of the cumulants do so approximately. Thus, although it is unlikely such cases will be encountered in inferential contexts, it is worth checking (by simulation) whether cumulants accurately replicate those of the target GQD. Furthermore, it may in some cases occur that the cumulants are indeed accurate whilst the density approximation fails. This can again be verified by simulation.

\bold{Warning [2]}:
The parameter \code{delt} is also used as the stepsize for solving a system of ordinary differential equations (ODEs) that govern the evolution of the cumulants of the diffusion. As such \code{delt} is required to be small for highly non-linear models in order to ensure sufficient accuracy.
}

\seealso{
See \code{\link{BiJGQD.mcmc}} and \code{\link{JGQD.density}}.
}



\examples{
\donttest{
#===============================================================================

 # For detailed notes and examples on how to use the BiJGQD.density() function, see
 # the following vignette:

 RShowDoc('Part_3_Bivariate_Diffusions',type='html','DiffusionRjgqd')

#===============================================================================
}
}
\keyword{ Transition density }
\keyword{ Moments }
