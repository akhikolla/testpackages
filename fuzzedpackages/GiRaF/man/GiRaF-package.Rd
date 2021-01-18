\name{GiRaF-package}
\alias{GiRaF-package}
\alias{GiRaF}
\docType{package}
\title{
  A toolbox for Gibbs Random Fields analysis
}
\description{
  \pkg{GiRaF} is a package for calculations on, and sampling from Gibbs (or discrete Markov) random fields.
}
\details{
\pkg{GiRaF} offers various tools for the analysis of Gibbs random fields and more precisely general homogeneous Potts model with possible anisotropy and potential on singletons (cliques composed of single vertex). \pkg{GiRaF} substantially lowers the barrier for practitioners aiming at analysing such Gibbs random fields. \pkg{GiRaF} contains exact methods for small lattices and several approximate methods for larger lattices that make the analysis easier for practitioners.

The \dQuote{GiRaF-introduction} vignette gives a detailled introduction on the package.

For a complete list of functions, use library(help = "GiRaF").
}
\author{
Julien Stoehr, Pierre Pudlo and Nial Friel.

Maintainer: Julien Stoehr <julien.stoehr@ucd.ie>
}
\references{
Friel, N. and Rue, H. (2007). Recursive computing and simulation-free inference for general factorizable models. \emph{Biometrika}, \bold{94(3)}:661--672.

Geman, S. and Geman, D. (1984). Stochastic Relaxation, Gibbs Distributions, and the Bayesian Restoration of Images. \emph{IEEE Transactions on Pattern Analysis and Machine Intellignence}, \bold{6(6)}:721-741.

Reeves, R. and Pettitt, A. N. (2004). Efficient recursions for general factorisable models. \emph{Biometrika}, \bold{91(3)}:751--757.

Swendsen, R. H. and Wang, J.-S. (1987). Nonuniversal critical dynamics in Monte Carlo simulations. \emph{Pysical Review Letters}, \bold{58(2)}:86-88.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{potts}
\keyword{mrf}
\keyword{gibbs}
\seealso{
The \dQuote{GiRaF-introduction} vignette
}
\examples{
# Dimension of the lattice
height <- 8
width <- 10

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

# Partition function for the default setting
NC.mrf(h = height, w = width, param = Beta)

# When specifying the number of colors and neighbors
NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta)

# When specifying an optional potential on sites
NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
       pot = potential)
       
# When specifying possible borders. The users will omit to mention all
# the non-existing borders
NC.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
       top = Top, left = Left, bottom = Bottom, right = Right, corner = Corner)
       
# Exact sampling for the default setting
exact.mrf(h = height, w = width, param = Beta, view = TRUE)

# When specifying the number of colors and neighbors
exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
          view = TRUE)
          
# When specifying an optional potential on sites
exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
       pot = potential, view = TRUE)
       
# When specifying possible borders. The users will omit to mention all
# the non-existing borders
exact.mrf(h = height, w = width, ncolors = K, nei = G, param = Beta, 
       top = Top, left = Left, bottom = Bottom, right = Right, corner = Corner, view = TRUE)
       
# Algorithm settings
n <- 200
method <- "Gibbs"

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
