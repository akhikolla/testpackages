\name{NC.mrf}
\alias{NC.mrf}
\docType{package}
\title{
Normalising constant of a Gibbs Random Field
}
\description{
Partition function of a general Potts model defined on a rectangular \eqn{h x w} lattice (\eqn{h \le w}) with either a first order or a second order dependency structure and a small number of rows (up to 25 for 2-state models).
}
\usage{
NC.mrf(h, w, param, ncolors = 2, nei = 4,  
       pot = NULL, top = NULL, left = NULL, 
       bottom = NULL, right = NULL, corner = NULL)
}
\arguments{
\item{h}{the number of rows of the rectangular lattice.}
\item{w}{the number of columns of the rectangular lattice.}
\item{param}{numeric entry setting the interaction parameter (edges parameter)}
\item{ncolors}{the number of states for the discrete random variables. By default, ncolors = 2.}
\item{nei}{the number of neighbors. The latter must be one of nei = 4 or nei = 8, which respectively correspond to a first order and a second order dependency structure. By default, nei = 4.}
\item{pot}{numeric entry setting homogeneous potential on singletons (vertices parameter). By default, pot = NULL}
\item{top, left, bottom, right, corner}{numeric entry setting constant borders for the lattice. By default, top = NULL, left = NULL, bottom = NULL, right = NULL, corner = NULL.}
}
\references{
Friel, N. and Rue, H. (2007). Recursive computing and simulation-free inference for general factorizable models. \emph{Biometrika}, \bold{94(3)}:661--672.

Reeves, R. and Pettitt, A. N. (2004). Efficient recursions for general factorisable models. \emph{Biometrika}, \bold{91(3)}:751--757.
}
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
}