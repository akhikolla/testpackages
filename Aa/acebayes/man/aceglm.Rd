\name{aceglm}
\alias{aceglm}
\alias{paceglm}

\title{
Approximate Coordinate Exchange (ACE) Algorithm for Generalised Linear Models
}
\encoding{UTF-8}
\description{
Functions implementing the approximate coordinate exchange (ACE) algorithm (Overstall & Woods, 2017) for finding Bayesian optimal experimental designs for generalised linear models (GLMs).
}
\usage{
aceglm(formula, start.d, family, prior, B, 
criterion = c("D", "A", "E", "SIG", "NSEL", "SIG-Norm", "NSEL-Norm"),
method = c("quadrature", "MC"), Q = 20, N1 = 20, N2 = 100, lower = -1, 
upper = 1, progress = FALSE, limits = NULL)

paceglm(formula, start.d, family, prior, B, 
criterion = c("D", "A", "E", "SIG", "NSEL", "SIG-Norm", "NSEL-Norm"),
method = c("quadrature", "MC"), Q = 20, N1 = 20, N2 = 100, lower = -1, 
upper = 1, limits = NULL, mc.cores = 1, n.assess = 20)
}

\arguments{
  \item{formula}{
An object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the column names of the argument \code{start.d}.
}
  \item{start.d}{
For \code{aceglm}, an \code{n} by \code{k} matrix, with column names used by the argument \code{formula}, specifying the initial design for the ACE algorithm.

For \code{paceglm}, a list with each element being an \code{n} by \code{k} matrix, with column names used by the argument \code{formula}, specifying the initial design for each repetition of the ACE algorithm.
}
	\item{family}{
A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
}
	\item{prior}{
An argument specifying the prior distribution.  	
	
For \code{method = "MC"}, a function with one argument: \code{B}; a scalar integer. This function should return a \code{B} by \code{p} matrix, with \code{p} the number of model parameters, containing a random sample from the prior distribution of the parameters. The value of \code{p} should correspond to the number of terms specified by the \code{formula} argument.

For \code{method = "quadrature"}, a list specifying a normal or uniform prior for the model parameters. For a normal prior distribution, the list should have named entries \code{mu} and \code{sigma2} specifying the prior mean and variance-covariance matrix. The prior mean may be specified as a scalar, which will then be replicated to form an vector of common prior means, or a vector of length \code{p}. The prior variance-covariance matrix may be specified as either a scalar common variance or a vector of length \code{p} of variances (for independent prior distributions) or as a \code{p} by \code{p} matrix. For a uniform prior distribution, the list should have a named entry \code{support}, a \code{2} by \code{p} matrix with each column giving the lower and upper limits of the support of the independent continuous uniform distribution for the corresponding parameter.
 
}
 \item{B}{
An optional argument for controlling the approximation to the expected utility. It should be a vector of length two. 

For \code{method = "MC"}, it specifies the size of the Monte Carlo samples, generated from the joint distribution of unknown quantities. The first sample size, \code{B[1]}, gives the sample size to use in the comparison procedures, and the second sample size, \code{B[2]}, gives the sample size to use for the evaluations of Monte Carlo integration that are used to fit the Gaussian process emulator. If left unspecified, the default value is \code{c(20000,1000)}.

For \code{method = "quadrature"}, it specifies the tuning parameters (numbers of radial abscissas and random rotations) for the implemented quadrature method; see \bold{Details} for more information. If left unspecified, the default value is \code{c(2, 8)}.
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
\bold{Shannon information gain} with normal-based Laplace approximation to marginal likelihood (\code{criterion = "SIG-Norm"});
\item
\bold{negative squared error loss} with importance sampling approximation to posterior mean (\code{criterion = "NSEL"});
\item
\bold{negative squared error loss} with normal-based approximation to posterior mean (\code{criterion = "NSEL-Norm"});
}
If left unspecified, the default is \code{"D"} denoting pseudo-Bayesian D-optimality. See \bold{Details} for more information. 
}
	\item{method}{
	An optional character argument specifying the method of approximating the expected utility function. Current choices are \code{method = "quadrature"} for a deterministic quadrature approximation and \code{method = "MC"} for a stochastic Monte Carlo approximation. The first of these choices is \bold{only} available when the argument \code{criterion = "A"}, \code{"D"} or \code{"E"}. The second choice is available for all possible values of the argument \code{criterion}. If left unspecified, the argument defaults to \code{"quadrature"} for \code{criterion = "A"}, \code{"D"} or \code{"E"} and to \code{"MC"} otherwise. See \bold{Details} for more information. 
}
 
  \item{Q}{
An integer specifying the number of evaluations of the approximate expected utility that are used to fit the Gaussian process emulator. The default value is \code{20}.
}
  \item{N1}{
An integer specifying the number of iterations of Phase I of the ACE algorithm (the coordinate exchange phase). 
The default value is \code{20}.
}
  \item{N2}{
An integer specifying the number of iterations of Phase II of the ACE algorithm (the point exchange phase).
The default value is \code{100}. 
}
  \item{lower}{
An argument specifying the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the lower limits of all coordinates of the design space. The default value is -1.
}
  \item{upper}{
An argument specifying the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the upper limits of all coordinates of the design space. The default value is 1.
}
  \item{progress}{
A logical argument indicating whether the iteration number and other information detailing the progress of the algorithm should be printed. The default value is \code{FALSE}.
}
  \item{limits}{
An argument specifying the grid over which to maximise the Gaussian process emulator for the expected utility function. It should be a function with three arguments: \code{i}, \code{j} and \code{d} which generates a one-dimensional grid for the \code{ij}th coordinate of the design when the current design is \code{d}. The default value is NULL which generates values uniformly on the interval \code{(lower[i,j],upper[i,j])} or \code{(lower,upper)} depending on whether the arguments \code{lower} and \code{upper} are matrices or scalars, respectively.
}
\item{mc.cores}{
The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one (the default), and parallelisation requires at least two cores. See \code{\link[parallel]{mclapply}} for more information and warnings for \code{mc.cores} > 1. 	
}
\item{n.assess}{
If \code{method = "MC"}, the approximate expected utility for the design from each repetition of the ACE algorithm will be calculated \code{n.assess} times. The terminal design returned will be the design with the largest mean approximate expected utility calculated over the \code{n.assess} approximations.}
}
\details{
The \code{aceglm} function implements the ACE algorithm to find designs for the class of generalised linear models (GLMs) for certain cases of utility function meaning the user does not have to write their own utility function.

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

Note that the utility functions for SIG and NSEL are currently only implemented for logistic regression, i.e. \code{family = binomial}, or Poisson regression, i.e. \code{family = poisson(link="log")}, whereas the utility functions for pseudo-Bayesian designs are implemented for generic GLM families.

Similar to all coordinate exchange algorithms, ACE should be repeated from different initial designs. The function 
\code{paceglm} will implement this where the initial designs are given by a list via the argument \code{start.d}. On the completion 
of the repetitions of ACE, \code{paceglm} will approximate the expected utility for all final designs and return the design (the terminal design) with the 
largest approximate expected utility.

For more details on the ACE algorithm, see Overstall & Woods (2017).
}

\value{
The function will return an object of class \code{"ace"} (for \code{aceglm}) or \code{"pace"} (for \code{paceglm})  which is a list with the following components:

\item{utility}{
The utility function resulting from the choice of arguments.
}
\item{start.d}{
The argument \code{start.d}.
}
\item{phase1.d}{
The design found from Phase I of the ACE algorithm.
}
\item{phase2.d}{
The design found from Phase II of the ACE algorithm.
}
\item{phase1.trace}{
A vector containing the evaluations of the approximate expected utility of the current design at each stage of Phase I of the ACE algorithm. This can be used to assess convergence.
}
\item{phase2.trace}{
A vector containing the evaluations of the approximate expected utility of the current design at each stage of Phase II of the ACE algorithm. This can be used to assess convergence.
}
\item{B}{
The argument \code{B}.	
}
\item{Q}{
The argument \code{Q}.
}
\item{N1}{
The argument \code{N1}.
}
\item{N2}{
The argument \code{N2}.
}
\item{glm}{
If the object is a result of a direct call to \code{aceglm} then this is \code{TRUE}.
}
\item{nlm}{
This will be \code{FALSE}.
}
\item{criterion}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{criterion}.
}
\item{method}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{method}.	
}
\item{prior}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{prior}.
}
\item{family}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{family}.
}
\item{formula}{
If the object is a result of a direct call to \code{acenlm} then this is the argument \code{formula}.	
}
\item{time}{
Computational time (in seconds) to run the ACE algorithm.
}
\item{binary}{
The argument \code{binary}. Will be \code{FALSE} for the utility functions currently implemented.
}
\item{d}{
The terminal design (\code{paceglm} only).
}
\item{eval}{
If \code{deterministic = "MC"}, a vector containing \code{n.assess} approximations to the expected utility for the terminal design (\code{paceglm} only).

If \code{deterministic = "quadrature"}, a scalar giving the approximate expected utility for the terminal design (\code{paceglm} only).
}
\item{final.d}{
A list of the same length as the argument \code{start.d}, where each element is the final design (i.e. \code{phase2.d}) for each repetition of the ACE algorithm (\code{paceglm} only).
}
\item{besti}{
A scalar indicating which repetition of the ACE algorithm resulted in the terminal design (\code{paceglm} only).
}
}

\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian experimental design: a review. \emph{Statistical
Science}, \bold{10}, 273-304.

Gotwalt, C. M., Jones, B. A.  & Steinberg, D. M. (2009). Fast computation of designs robust to parameter uncertainty for nonlinear settings. \emph{Technometrics}, \bold{51}, 88-95.

Meyer, R. & Nachtsheim, C. (1995). The coordinate exchange algorithm for constructing exact
optimal experimental designs. \emph{Technometrics}, \bold{37}, 60-69.

Monahan, J. and Genz, A. (1997). Spherical-radial integration rules for
Bayesian computation,” Journal of the American Statistical Association, 92, 664–674.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.

Overstall, A.M., McGree, J.M. & Drovandi, C.C. (2018). An approach for finding fully Bayesian optimal designs using normal-based approximations to loss functions. \emph{Statistics and Computing}, \bold{28}(2), 343-358.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\note{
These are wrapper functions for \code{\link{ace}} and \code{\link{pace}}.
}

\seealso{
\code{\link{ace}}, \code{\link{acenlm}}, \code{\link{pace}}, \code{\link{pacenlm}}.
}
\examples{
## This example uses aceglm to find a Bayesian D-optimal design for a 
## first order logistic regression model with 6 runs 4 factors. The priors are 
## those used by Overstall & Woods (2017), with each of the five
## parameters having a uniform prior. The design space for each coordinate is [-1, 1].

set.seed(1)
## Set seed for reproducibility.

n<-6
## Specify the sample size (number of runs).

start.d<-matrix(2 * randomLHS(n = n,k = 4) - 1,nrow = n,ncol = 4,
dimnames = list(as.character(1:n), c("x1", "x2", "x3", "x4")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.

low<-c(-3, 4, 5, -6, -2.5)
upp<-c(3, 10, 11, 0, 3.5)
## Lower and upper limits of the uniform prior distributions.

prior<-function(B){
t(t(6*matrix(runif(n = 5 * B),ncol = 5)) + low)}
## Create a function which specifies the prior. This function will return a 
## B by 5 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.

example1<-aceglm(formula=~x1+x2+x3+x4, start.d = start.d, family = binomial, 
prior = prior, method = "MC", N1 = 1, N2 = 0, B = c(1000, 1000))
## Call the aceglm function which implements the ACE algorithm requesting 
## only one iteration of Phase I and zero iterations of Phase II. The Monte
## Carlo sample size for the comparison procedure (B[1]) is set to 100.

example1
## Print out a short summary.

#Generalised Linear Model 
#Criterion = Bayesian D-optimality 
#Formula: ~x1 + x2 + x3 + x4
#
#Family: binomial 
#Link function: logit 
#
#Method:  MC 
#
#B:  1000 1000 
#
#Number of runs = 6
#
#Number of factors = 4
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 0
#
#Computer time = 00:00:01

example1$phase2.d
## Look at the final design.

#          x1          x2          x3         x4
#1 -0.4735783  0.12870470 -0.75064318  1.0000000
#2 -0.7546841  0.78864527  0.58689270  0.2946728
#3 -0.7463834  0.33548985 -0.93497463 -0.9573198
#4  0.4446617 -0.29735212  0.74040030  0.2182800
#5  0.8459424 -0.41734194 -0.07235575 -0.4823212
#6  0.6731941  0.05742842  1.00000000 -0.1742566

prior2 <- list(support = rbind(low, upp))
## A list specifying the parameters of the uniform prior distribution

example2<-aceglm(formula = ~ x1 +x2 + x3 + x4, start.d = start.d, family = binomial, 
prior = prior2, N1 = 1, N2 = 0)
## Call the aceglm function with the default method of "quadrature"

example2$phase2.d
## Final design

#          x1          x2          x3         x4
#1 -0.4647271  0.07880018 -0.94648750  1.0000000
#2 -0.7102715  0.79827332  0.59848578  0.5564422
#3 -0.7645090  0.39778176 -0.74342036 -1.0000000
#4  0.4514632 -0.33687477  0.55066110  0.3994593
#5  0.7913559 -0.41856994  0.01321035 -0.8848135
#6  0.6337306  0.11578522  1.00000000  1.0000000

mean(example1$utility(d = example1$phase2.d, B = 20000))
#[1] -11.61105
mean(example2$utility(d = example2$phase2.d, B = 20000))
#[1] -11.19737
## Compare the two designs using the Monte Carlo approximation
}
