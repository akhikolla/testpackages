
<!-- README.md is generated from README.Rmd. Please edit that file -->
Rdtq
====

Rdtq implements Density Tracking by Quadrature (DTQ) algorithms for stochastic differential equations. Suppose you have a stochastic differential equation of the form

dX(t) = f(X(t)) dt + g(X(t)) dW(t)

where W(t) is standard Brownian motion, and where f and g are user-specified drift and diffusion functions.

Let p(x,t) denote the probability density function of the random variable X(t). Then Rdtq will calculate a numerical approximation to p(x,T) for a fixed time T &gt; 0.

There are three main reasons to use Rdtq to compute the density function of X(t):

1.  *Rdtq is simultaneously easy to use and flexible*. The only functions the user has to write are those that return the drift f(x) and the diffusion g(x). In particular, the user does not have to specify the derivatives of these functions. Beyond these two functions, the user can specify the initial condition either in the form of a single point, X\_0 = C, or in the form of a density function, p(x,0) at each grid point x. The user also has total control over the temporal and spatial grid spacings.

2.  *Rdtq is provably convergent*, in the limit where the temporal and spatial grid spacings vanish. See [H. S. Bhat and R. W. M. A. Madushani, "Density Tracking by Quadrature for Stochastic Differential Equations," arXiv:1610.09572](http://bit.ly/2fbNsp5).

3.  *Rdtq is fast*. The package includes both C++ and sparse matrix implementations of the DTQ algorithm. In our tests at the finest grid resolutions, the DTQ method is approximately 100 times faster than a competing method that consists of numerically solving the Fokker-Planck/Kolmogorov partial differential equation.

Installation
------------

You can install Rdtq from github with:

``` r
# install.packages("devtools")
devtools::install_github("hbhat4000/Rdtq")
```

Example
-------

Suppose you have the stochastic differential equation dX(t) = -X(t) dt + dW(t) and you would like to find the probability density function of X(t) at T=1. Here is how to do that using Rdtq.

``` r
require(Rdtq)
#> Loading required package: Rdtq
# We use the drift function f(x) = -x and diffusion function g(x) = 1.
mydrift = function(x) { -x }
mydiff = function(x) { rep(1,length(x)) }
# We use the sparse matrix implementation of DTQ and solve for the density
# function at final time fT=1.  The initial condition is X(0)=0,
# the temporal step size is h=0.1, the grid spacing is k=0.01, and 
# the spatial grid extends from -250*0.01 to 250*0.01.
test = rdtq(h=0.1,k=0.01,bigm=250,init=0,fT=1,
            drift=mydrift,diffusion=mydiff,method="sparse")
#> [1] "Calling dtq with grid specified via k and bigm."
#> [1] "Using sparse method."
plot(test$xvec,test$pdf,type='l')
```

![](README-example-1.png)
