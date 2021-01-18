\name{sampler.mrf}
\alias{sampler.mrf}
\docType{package}
\title{
MCMC samplers for Gibbs Random Fields
}
\description{
sampler.mrf gives approximate sample from the likelihood of a general Potts model defined on a rectangular \eqn{h x w} lattice (\eqn{h \le w}) with either a first order or a second order dependency structure. Available options are the Gibbs sampler (Geman and Geman (1984)) and the Swendsen-Wang algorithm (Swendsen and Wang (1987)).
}
\usage{
sampler.mrf(iter, sampler = "Gibbs" , h, w, 
          param, ncolors = 2, nei = 4, pot = NULL, 
          top = NULL, left = NULL, bottom = NULL, right = NULL, 
          corner = NULL, initialise = TRUE, random = TRUE, view = FALSE)
}
\arguments{
\item{iter}{Number of iterations of the algorithm.}
\item{sampler}{The method to be used. The latter must be one of "Gibbs" or "SW" corresponding respectively to the Gibbs sampler and the Swendsen-Wang algorithm.}
\item{h}{the number of rows of the rectangular lattice.}
\item{w}{the number of columns of the rectangular lattice.}
\item{param}{numeric entry setting the interaction parameter (edges parameter)}
\item{ncolors}{the number of states for the discrete random variables. By default, ncolors = 2.}
\item{nei}{the number of neighbors. The latter must be one of nei = 4 or nei = 8, which respectively correspond to a first order and a second order dependency structure. By default, nei = 4.}
\item{pot}{numeric entry setting homogeneous potential on singletons (vertices parameter). By default, pot = NULL}
\item{top, left, bottom, right, corner}{numeric entry setting constant borders for the lattice. By default, top = NULL, left = NULL, bottom = NULL, right = NULL, corner = NULL.}
\item{initialise}{Logical value indicating whether initial guess should be randomly drawn.}
\item{random}{Logical value indicating whether the sites should be updated sequentially or randomdly. Used only with the "Gibbs" option.}
\item{view}{Logical value indicating whether the draw should be printed. Do not display the optional borders.}
}
\references{
Geman, S. and Geman, D. (1984). Stochastic Relaxation, Gibbs Distributions, and the Bayesian Restoration of Images. \emph{IEEE Transactions on Pattern Analysis and Machine Intellignence}, \bold{6(6)}:721-741.

Swendsen, R. H. and Wang, J.-S. (1987). Nonuniversal critical dynamics in Monte Carlo simulations. \emph{Pysical Review Letters}, \bold{58(2)}:86-88.
}
\seealso{
The \dQuote{GiRaF-introduction} vignette
}
\examples{
# Algorithm settings
n <- 200
method <- "Gibbs"

# Dimension of the lattice
height <- width <- 100

# Interaction parameter
Beta <- 0.6 # Isotropic configuration
# Beta <- c(0.6, 0.6) # Anisotropic configuration when nei = 4
# Beta <- c(0.6, 0.6, 0.6, 0.6) # Anisotropic configuration when nei = 8

# Number of colors
K <- 2 
# Number of neighbors 
G <- 4

# Optional potential on sites
potential <- runif(K,-1,1)
# Optional borders. 
Top <- Bottom <- sample(0:(K-1), width, replace = TRUE)
Left <- Right <- sample(0:(K-1), height, replace = TRUE)
Corner <- sample(0:(K-1), 4, replace = TRUE)

# Sampling method for the default setting
sampler.mrf(iter = n, sampler = method, h = height, w = width, 
            param = Beta, view = TRUE)
            
# Sampling using an existing configuration as starting point
sampler.mrf(iter = n, sampler = method, h = height, w = width, 
            ncolors = K, nei = G, param = Beta, 
            initialise = FALSE, view = TRUE)
            
# Specifying optional arguments. The users may omit to mention all
# the non-existing borders
sampler.mrf(iter = n, sampler = method, h = height, w = width, 
            ncolors = K, nei = G, param = Beta,
            pot = potential, top = Top, left = Left, bottom = Bottom, 
            right = Right, corner = Corner, view = TRUE)
            
# Gibbs sampler with sequential updates of the sites. 
sampler.mrf(iter = n, sampler = "Gibbs", h = height, w = width, 
            ncolors = K, nei = G, param = Beta,
            random = FALSE, view = TRUE)
}