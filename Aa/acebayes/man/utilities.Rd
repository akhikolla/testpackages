\name{utilities}
\alias{utilitynlm}
\alias{utilityglm}

\title{
Approximate expected utility function for generalised linear models and non-linear regression models
}
\encoding{UTF-8}
\description{
Generates an approximate utility function for generalised linear models and non-linear regression models.
}
\usage{
utilityglm(formula, family, prior, 
criterion = c("D", "A", "E", "SIG", "NSEL", "SIG-Norm", "NSEL-Norm"),
method = c("quadrature", "MC"), nrq)

utilitynlm(formula, prior, desvars, criterion =  c("D", "A", "E", "SIG", "NSEL"),
method = c("quadrature", "MC"), nrq)
}

\arguments{
  \item{formula}{
An argument providing a symbolic description of the model.

For \code{utilityglm}, an object of class \code{"formula"}: a symbolic description of the model. 

For \code{utilitynlm}, an object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the argument \code{prior}
}
	\item{family}{
For \code{utilityglm}, a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
}
	\item{prior}{
An argument specifying the prior distribution.  	
	
For \code{method = "MC"}, a function with one argument: \code{B}; a scalar integer. This function should return a \code{B} by \code{p} matrix (\code{p+1} for \code{criterion = "SIG"} or \code{criterion = "NSEL"}), with \code{p} the number of model parameters, containing a random sample from the prior distribution of the parameters. The value of \code{p} should correspond to the number of terms specified by the \code{formula} argument. For \code{utilitynlm}, the column names must match the names of parameters in the \code{formula} argument. For \code{utilitynlm}, if \code{criterion="SIG"}, \code{criterion="NSEL"}, \code{criterion="SIG-Norm"} or \code{criterion="NSEL-Norm"} then an extra column called \code{sig2} should be included with a sample from the error variance.

For \code{method = "quadrature"}, a list specifying a normal or uniform prior for the model parameters. For a normal prior distribution, the list should have named entries \code{mu} and \code{sigma2} specifying the prior mean and variance-covariance matrix. The prior mean may be specified as a scalar, which will then be replicated to form an vector of common prior means, or a vector of length \code{p}. The prior variance-covariance matrix may be specified as either a scalar common variance or a vector of length \code{p} of variances (for independent prior distributions) or as a \code{p} by \code{p} matrix. For \code{utilitynlm}, the names attribute of \code{mu} must match the names of the parameters in the \code{formula} argument. For a uniform prior distribution, the list should have a named entry \code{support}, a \code{2} by \code{p} matrix with each column giving the lower and upper limits of the support of the independent continuous uniform distribution for the corresponding parameter. For \code{utilitynlm}, the column names of \code{support} must match the names of parameters in the \code{formula} argument.	
}

\item{desvars}{
For \code{utilitynlm}, a character vector listing the design variables that appear in the argument \code{formula}.
}
	\item{criterion}{
An optional character argument specifying the utility function. There are currently seven utility functions implemented as follows:
\enumerate{
\item
\bold{pseudo-Bayesian D-optimality} (\code{criterion = "D"});
\item
\bold{pseudo-Bayesian A-optimality} (\code{criterion = "A"});
\item
\bold{pseudo-Bayesian E-optimality} (\code{criterion = "E"}).
\item
\bold{Shannon information gain} with Monte Carlo (MC) approximation to marginal likelihood (\code{criterion = "SIG"});
\item
\bold{Shannon information gain} with normal-based Laplace approximation to marginal likelihood (\code{criterion = "SIG-Norm"}, only for \code{utilityglm}));
\item
\bold{negative squared error loss} with importance sampling approximation to posterior mean (\code{criterion = "NSEL"});
\item
\bold{negative squared error loss} with normal-based approximation to posterior mean (\code{criterion = "NSEL-Norm"}, only for \code{utilityglm})) ;
}
The default value is \code{"D"} denoting pseudo-Bayesian D-optimality. See \bold{Details} for more information. 
}
\item{method}{
	An optional character argument specifying the method of approximating the expected utility function. Current choices are \code{method = "quadrature"} for a deterministic quadrature approximation and \code{method = "MC"} for a stochastic Monte Carlo approximation. The first of these choices is \bold{only} available when the argument \code{criterion = "A"}, \code{"D"} or \code{"E"}. The second choice is available for all possible values of the argument \code{criterion}. If left unspecified, the argument defaults to \code{"quadrature"} for \code{criterion = "A"}, \code{"D"} or \code{"E"} and to \code{"MC"} otherwise. See \bold{Details} for more information. 			
}
\item{nrq}{
For \code{method = "quadrature"}, a vector of length two specifying the number of radial abscissas (\code{nrq[1]}) and quasi-random rotations (\code{nrq[2]}) required for the implemented quadrature scheme; see \bold{Details} for more information. If left unspecified, the default value is \code{c(2, 8)}.	
}
}
\details{
Two utility functions are implemented.
\enumerate{
\item
\bold{Shannon information gain (SIG)} 

The utility function is
\deqn{u^{SIG}(d) = \pi(\theta|y,d) - \pi(\theta),}{
U^SIG(d) = f(\theta|y,d) - f(\theta),}
where \eqn{\pi(\theta|y,d)}{f(\theta|y,d)} and \eqn{\pi(\theta)}{f(\theta)} denote the posterior and prior densities of the parameters \eqn{\theta}{\theta}, respectively.
\item
\bold{Negative squared error loss (NSEL)}

The utility function is
\deqn{u^{NSEL}(d) = - \left(\theta - E(\theta |y,d)\right)^T \left(\theta - E(\theta |y,d)\right),}{
u^NSEL(d) = - (\theta - E(\theta |y,d))^T(\theta - E(\theta |y,d)),}
where \eqn{E(\theta | y,d)}{E(\theta | y,d)} denotes the posterior mean of \eqn{\theta}{\theta}.
}

In both cases the utility function is not available in closed form due to the analytical intractability of either the posterior distribution (for SIG) or the posterior mean (for NSEL). The \code{acebayes} package implements two approximations to both utility functions. If \code{criterion = "SIG"} or \code{criterion = "NSEL"} then sampling-based Monte Carlo or importance sampling approximations will be employed. This was the original approach used by Overstall & Woods (2017). If \code{criterion = "SIG-Norm"} or \code{criterion = "NSEL-Norm"} then approximations based on approximate normality of the posterior (Overstall et al., 2017) will be used.

The normal approximation to the posterior can be taken further leading to the approximation by some scalar function of the Fisher information matrix, \eqn{\mathcal{I} (\theta;d)}{I (\theta;d)}, which only depends on \eqn{\theta}{\theta} (Chaloner & Verdinelli, 1995). In the case of SIG, the approximate utility is given by
\deqn{u^{D}(d) = \log \vert \mathcal{I} (\theta;d) \vert,}{
u^D(d) = log | I(\theta;d)|,}
and the resulting design is typically called pseudo-Bayesian D-optimal. For NSEL, the approximate utility is given by
\deqn{u^A(d) = - \mathrm{tr} \left\{ \mathcal{I} (\theta;d)^{-1} \right\} }{
u^A(d) = - tr (I(\theta;d)^(-1))}
with the resulting design termed pseudo-Bayesian A-optimal. These designs are often used under the frequentist approach to optimal experimental design and so to complete the usual set, the following utility for finding a pseudo-Bayesian E-optimal design is also implemented:
\deqn{u^E(d) = \mathrm{min} \mbox{ } e\left(\mathcal{I} (\theta;d) \right),}{
U^E(d) = min(e(I(\theta;d))),}
where \eqn{e()}{e()} denotes the function that calculates the eigenvalues of its argument.

The expected utilities can be approximated using Monte Carlo methods (\code{method = "MC"} for all criteria) or using a deterministic quadrature method (\code{method = "quadrature"}, implemented for the D, A and E criteria). The former approach approximates the expected utility via sampling from the prior. The latter approach uses a radial-spherical integration rule (Monahan and Genz, 1997) and \code{B[1]} specifies the number, \eqn{n_r}{n_r}, of radial abscissas and \code{B[2]} specifies the number,  \eqn{n_q}{n_q}, of random rotations. Larger values of  \eqn{n_r}{n_r} will produce more accurate, but also more computationally expensive, approximations. See Gotwalt et al. (2009) for further details.   

For \code{utilityglm}, note that the utility functions for SIG and NSEL are currently only implemented for logistic regression, i.e. \code{family = binomial}, or Poisson regression, i.e. \code{family = poisson(link = "log")}, whereas the utility functions for pseudo-Bayesian designs are implemented for generic GLM families.

For more details on the ACE algorithm, see Overstall & Woods (2017).
}

\value{
The function will return a list with the following components:
\item{utility}{
The utility function resulting from the choice of arguments.
}
}

\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian experimental design: a review. \emph{Statistical
Science}, \bold{10}, 273-304.

Gotwalt, C. M., Jones, B. A.  & Steinberg, D. M. (2009). Fast computation of designs robust to parameter uncertainty for nonlinear settings. \emph{Technometrics}, \bold{51}, 88-95.

Monahan, J. and Genz, A. (1997). Spherical-radial integration rules for
Bayesian computation,” Journal of the American Statistical Association, 92, 664-674.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\seealso{
\code{\link{aceglm}}, \code{\link{acenlm}}, \code{\link{paceglm}}, \code{\link{pacenlm}}.
}
\examples{
## 1. This example uses utilityglm to generate the pseudo-Bayesian D-optimality
## approximate expected utility function using a Monte Carlo approximation.

low<-c(-3, 4, 5, -6, -2.5)
upp<-c(3, 10, 11, 0, 3.5)
## Lower and upper limits of the uniform prior distributions.

prior<-function(B){
t(t(6*matrix(runif(n=5*B),ncol=5))+low)}
## Create a function which specifies the prior. This function will return a 
## B by 5 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.

ex <- utilityglm(formula = ~x1+x2+x3+x4, family = binomial, prior = prior, method = "MC")

set.seed(1)
## Set seed for reproducibility.

n<-6
## Specify the sample size (number of runs).

start.d<-matrix( 2 * randomLHS(n = n,k = 4) - 1,nrow = n,ncol = 4,
dimnames = list(as.character(1:n),c("x1", "x2", "x3", "x4")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.

ex$utility(d = start.d, B = 10)
## Evaluate resulting approximate utility. Should get:

#[1] -13.98143 -17.07772 -19.88988 -22.40720 -15.27411 -15.02717 -16.17253 -18.66600 -13.75118
#[10] -21.83820

## 2. This example uses utilitynlm to generate the psuedo-Bayesian A-optimality expected utility
## function using a quadrature approximation

low<-c(0.01884, 0.298, 21.8)
upp<-c(0.09884, 8.298, 21.8)
## Lower and upper limits of the uniform prior distributions. Note that the prior
## for the third element is a point mass.

prior2 <- list(support = cbind(rbind(low, upp)))
colnames(prior2$support) <- c("a", "b", "c")
## Specify a uniform prior with ranges given by low and upp

ex2 <- utilitynlm(formula = ~ c * (exp( - a * t) - exp( - b *t)), prior = prior2, 
                            desvars = "t")
                            
n <- 6
start.d <- matrix(24 * randomLHS(n = n, k = 1), nrow = n)
colnames(start.d) <- "t"
ex2$utility(d = start.d) 
## -13.17817    

}
