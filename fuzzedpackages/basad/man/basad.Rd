\name{basad}
\alias{basad}
\title{Bayesian variable selection with shrinking and diffusing priors}
\description{
    This function performs the Bayesian variable selection procedure with shrinking and diffusing priors via Gibbs sampling. Three different prior options placed on the coefficients are provided: Gaussian, Student's t, Laplace. The posterior estimates of coefficients are returned and the final model is selected either by using the "BIC" criterion or the median probability model.
}
\usage{
basad( x = NULL, y = NULL, K = -1, df = 5, nburn = 1000, niter = 1000,
        alternative = FALSE, verbose = FALSE, nsplit = 20, tau0 = NULL, tau1 = NULL,
        prior.dist = "Gauss", select.cri = "median", BIC.maxsize = 20)
}
\arguments{
\item{x}{The matrix or data frame of covariates. }
\item{y}{The response variables. }
\item{K}{An initial value for the numbers of active covariates in the model. This value is related to the prior probability that a covariate is nonzero. If \code{K} is not specified greater than 3, this prior probability will be estimated by a Beta prior using Gibbs sampling (see details below).  }
\item{df}{The degrees of freedom of t prior when \code{prior.dist == "t"}. }
\item{nburn}{The number of iterations for burn-in. }
\item{niter}{The number of iterations for estimation. }
\item{alternative}{If \code{TRUE}, an alternative sampling scheme from Bhattacharya will be used which can accelerate the speed of the algorithm for very large p. However, when using block updating (by setting \code{nsplit} to be greater than 1) this alternative sampling will not be invoked.  }
\item{verbose}{If \code{TRUE}, verbose output is sent to the terminal. }
\item{nsplit}{Numbers of splits for the block updating scheme. }
\item{tau0}{ The scale of the prior distribution for inactive coefficients (see details below). }
\item{tau1}{ The scale of the prior distribution for active coefficients (see details below). }
\item{prior.dist}{Choice of the base distribution for spike and slab priors. If \code{prior.dist="t"}, the algorithm will place Student's t prior for regression coefficients. If \code{prior.dist="Laplace"}, the algorithm will place Laplace prior. Otherwise, it will place the default Gaussian priors. }
\item{select.cri}{Model selection criteria. If \code{select.cri="median"}, the algorithm will use the median probability model to select the active variables. If \code{ select.cri="BIC"}, the algorithm will use the BIC criteria to select the active variables. }
\item{BIC.maxsize}{ The amount of the variables that are chosen to apply BIC criteria based on the ranking of their marginal posterior probabilities. If the input sample size is less than the default value 20, all variables will be considered when applying BIC. }
}


\value{

An object of class \code{basad} with the following components:
\item{all.var}{ Summary object for all the variables. }
\item{select.var}{ Summary object for the selected variables. }
\item{beta.names}{ Variable names for the coefficients. }
\item{verbose}{ Verbose details (used for printing). }
\item{posteriorZ}{A vector of the marginal posterior probabilities for the latent vector Z. }
\item{model.index}{A vector containing the indices of selected variables. }
\item{modelZ}{A binary vector Z indicating whether the coefficient is true in the selected model. }
\item{est.B}{ Estimated coefficient values from the posterior distribution through Gibbs sampling. }
\item{allB}{ A matrix of all sampled coefficient values along the entire chain. Each row represents the sampled values under each iteration.}
\item{allZ}{ A matrix of all sampled posterior probabilities for the latent variable Z along the entire chain. Each row represents the sampled values under each iteration. }
\item{x}{ Standardized x-matrix. }
\item{y}{ Standardized y vector. }
}
\details{


In the package, the regression coefficients have following hierarchical structure:
\deqn{ \beta | ( Z = 0, \sigma^2 )  = N(0, \tau_0^2 \sigma^2 ),  \beta | ( Z = 1, \sigma^2 )  = N(0,  \tau_1^2 \sigma^2 )}

where the latent variable \eqn{Z_i} of value 0 or 1 indicates whether \eqn{i}th variable is in the slab and spike part of the prior. The package provides different prior choices for the coefficients: Gaussian, Student's t, Laplace. Through setting the parameter \code{prior.dist}, the coefficients will have the corresponding prior densities as follows:

1. The Gaussian priors case:
\deqn{ \beta | ( Z = k, \sigma^2 ) =  \frac{1}{ \sqrt{ 2 \pi \tau^2_k \sigma^2} }  e^{ - \frac{   \beta^2  }{ 2 \tau^2_k \sigma^2  }    }  }

2. The Student's t prior case:
\deqn{ \beta | ( Z = k, \sigma^2 ) = \frac{  \Gamma ( \frac{\nu + 1 }{2} )  }{ \Gamma ( \frac{\nu}{2} ) \sqrt{ \pi \nu } \tau_k \sigma } \Bigg(  1 + \frac{1}{\nu}  \Big( \frac{ \beta^2 }{ \tau^2_k \sigma^2 } \Big)   \Bigg)^{  -\frac{ \nu + 1}{ 2 } } }

Where \eqn{\nu} is the degrees of freedom

3. The Laplace prior case:
\deqn{  \beta | ( Z = k, \sigma^2 )  = \frac{1}{2 \tau_k^2 \sigma^2 } \exp \Big(  -\frac{  |\beta|  }{  \tau_{k}^2 \sigma^2  }  \Big) }

The \eqn{\tau_k} is the scale for the prior distribution. If user did not set a specific value, the prior scales are specified as follows:
\deqn{
\tau_{0}^2 = \frac{1}{n} a_\tau , \,  \tau_{1}^2 = \max\Big( 100 \tau_{0}^2, \frac{ \tau_{0}  p_n }{ (1 - p_n) \rho  } \Big),
}

where \eqn{ \rho} is the prior density evaluated at \eqn{ f_p( b_\tau \times \log{(p_n + 1)})}, \eqn{ f_p } is the density function for the corresponding prior distribution. The parameter \eqn{a} and \eqn{b} are \eqn{a_\tau = 1} and \eqn{b_\tau = 2.4} by default.


The prior probability \eqn{q_n = P(Z_i = 1)} that a covariate is nonzero can be specified by value K. The K represents a prior belief of the upper bound of the true covariates in the model. When user specifies a value of K greater than 3, setting \eqn{q_n = c/p_n}, through the calculation(see details in Naveen (2014)):
\deqn{ \Phi( (K - c)/ \sqrt{c} ) = 1 - \alpha }
The prior probability on the models with sizes greater than K will be \eqn{\alpha}, and this \eqn{\alpha} is set to 0.1 in the package.


}

\references{
Narisetty, N. N., & He, X. (2014). Bayesian variable selection with shrinking and diffusing priors. \emph{The Annals of Statistics}, 42(2), 789-817.

Bhattacharya, A., Chakraborty, A., & Mallick, B. K. (2016). Fast sampling with Gaussian scale mixture priors in high-dimensional regression. \emph{Biometrika}, 4(103), 985-991.

Barbieri, M. M., & Berger, J. O. (2004). Optimal predictive model selection. \emph{The Annals of Statistics}, 32(3), 870-897.
}
\examples{
\dontshow{
#-----------------------------------------------------------
#Generate Data: The simulated low dimensional data
#-----------------------------------------------------------

n =  100; p = 99; nz = 5

rho1=0.25;rho2=0.25;rho3=0.25  ### correlations
Bc = c( 0,seq(0.6,3,length.out=nz), array(0, p-nz))

covr1=(1- rho1)*diag(nz) +  array(rho1,c(nz,nz))
covr3=(1- rho3)*diag(p-nz) +  array(rho3,c(p-nz,p-nz))
covr2=array(rho2,c(nz,p-nz))
covr=rbind( cbind(covr1,covr2), cbind(t(covr2),covr3) )

covE = eigen(covr)
covsq = covE$vectors \%*\% diag( sqrt(covE$values) ) \%*\% t(covE$vectors)

Xs = matrix( rnorm(n*p), nrow = n); Xn = covsq \%*\% t(Xs)
X = cbind(array(1, n), t(Xn))
Y = X \%*\% Bc + rnorm(n); X <- X[,2:ncol(X)]


#-----------------------------------------------------------
#Example 1: Run the default setting of the Guassian priors
#-----------------------------------------------------------

obj <- basad( x = X, y = Y)
print( obj )
}
\donttest{

#-----------------------------------------------------------
#Generate Data: The simulated high dimensional data
#-----------------------------------------------------------

n =  100; p = 499; nz = 5

rho1=0.25;rho2=0.25;rho3=0.25  ### correlations
Bc = c( 0,seq(0.6,3,length.out=nz), array(0, p-nz))

covr1=(1- rho1)*diag(nz) +  array(rho1,c(nz,nz))
covr3=(1- rho3)*diag(p-nz) +  array(rho3,c(p-nz,p-nz))
covr2=array(rho2,c(nz,p-nz))
covr=rbind( cbind(covr1,covr2), cbind(t(covr2),covr3) )

covE = eigen(covr)
covsq = covE$vectors \%*\% diag( sqrt(covE$values) ) \%*\% t(covE$vectors)

Xs = matrix( rnorm(n*p), nrow = n); Xn = covsq \%*\% t(Xs)
X = cbind(array(1, n), t(Xn))
Y = X \%*\% Bc + rnorm(n); X <- X[,2:ncol(X)]


#-----------------------------------------------------------
#Example 1: Run the default setting of the Guassian priors
#-----------------------------------------------------------

obj <- basad( x = X, y = Y)
print( obj )

#-----------------------------------------------------------
#Example 2: Use different priors and slection criteria
#-----------------------------------------------------------

obj <- basad( x = X, y = Y, prior.dist = "t", select.cri = "BIC")
print( obj )
}



}
\author{
Qingyan Xiang (\email{qxiang@illinois.edu})

Naveen Narisetty (\email{naveen@illinois.edu})
}
\keyword{regression}

