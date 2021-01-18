\name{acenlm}
\alias{acenlm}
\alias{pacenlm}

\title{
Approximate Coordinate Exchange (ACE) Algorithm for Non-Linear Models
}
\encoding{UTF-8}
\description{
Functions implementing the approximate coordinate exchange algorithm (Overstall & Woods, 2017) for finding optimal Bayesian designs for non-linear regression models.
}
\usage{
acenlm(formula, start.d, prior, B, criterion = c("D", "A", "E", "SIG", "NSEL"), 
method = c("quadrature", "MC"),  Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
progress = FALSE, limits = NULL)

pacenlm(formula, start.d, prior, B, criterion = c("D", "A", "E", "SIG", "NSEL"), 
method = c("quadrature", "MC"),  Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
limits = NULL, mc.cores = 1, n.assess = 20)
}

\arguments{
  \item{formula}{
An object of class \code{"formula"}: a symbolic description of the model. The terms should correspond to the column names of the argument \code{start.d} and the argument \code{prior}.
}
  \item{start.d}{
For \code{aceglm}, an \code{n} by \code{k} matrix, with column names used by the argument \code{formula}, specifying the initial design for the ACE algorithm.

For \code{paceglm}, a list with each element, an \code{n} by \code{k} matrix, with column names used by the argument \code{formula}, specifying the initial design for each repetition of the ACE algorithm.
}
	\item{prior}{
An argument specifying the prior distribution.  	
	
For \code{method = "MC"}, a function with one argument: \code{B}; a scalar integer. This function should return a \code{B} by \code{p} matrix (or \code{p+1} for \code{criterion = "SIG"} or \code{criterion = "NSEL"}), with \code{p} the number of model parameters, containing a random sample from the prior distribution of the parameters. The value of \code{p} should correspond to the number of terms specified by the \code{formula} argument. The column names must match the names of parameters in the \code{formula} argument. For \code{criterion="SIG"} or \code{criterion="NSEL"}, an extra column (named \code{sig2}) should contain a sample from the prior distribution of the error variance.

For \code{method = "quadrature"}, a list specifying a normal or uniform prior for the model parameters. For a normal prior distribution, the list should have named entries \code{mu} and \code{sigma2} specifying the prior mean and variance-covariance matrix. The prior mean may be specified as a scalar, which will then be replicated to form an vector of common prior means, or a vector of length \code{p}. The prior variance-covariance matrix may be specified as either a scalar common variance or a vector of length \code{p} of variances (for independent prior distributions) or as a \code{p} by \code{p} matrix. The names attribute of \code{mu} must match the names of the parameters in the \code{formula} argument. For a uniform prior distribution, the list should have a named entry \code{support}, a \code{2} by \code{p} matrix with each column giving the lower and upper limits of the support of the independent continuous uniform distribution for the corresponding parameter. The column names of \code{support} must match the names of parameters in the \code{formula} argument.
}
 \item{B}{
An optional argument for controlling the approximation to the expected utility. It should be a vector of length two. 

For \code{criterion = "MC"}, it specifies the size of the Monte Carlo samples, generated from the joint distribution of unknown quantities. The first sample size, \code{B[1]}, gives the sample size to use in the comparison procedures, and the second sample size, \code{B[2]}, gives the sample size to use for the evaluations of Monte Carlo integration that are used to fit the Gaussian process emulator. If left unspecified, the default value is \code{c(20000,1000)}.

For \code{criterion = "quadrature"}, it specifies the tuning parameters (numbers of radial abscissas and random rotations) for the implemented quadrature method; see \bold{Details} for more information. If left unspecified, the default value is \code{c(2, 8)}.
}

	\item{criterion}{
An optional character argument specifying the utility function. There are currently five utility functions implemented consisting of 
\enumerate{
\item
\bold{pseudo-Bayesian D-optimality} (\code{criterion = "D"}); 
\item
\bold{pseudo-Bayesian A-optimality} (\code{criterion = "A"}); 
\item
\bold{pseudo-Bayesian E-optimality} (\code{criterion = "E"}); 
\item
\bold{Shannon information gain} (\code{criterion = "SIG"}); 
\item
\bold{negative squared error loss} (\code{criterion = "NSEL"}). 
}
The default value is \code{"D"} denoting pseudo-Bayesian D-optimality. See \bold{Details} for more information. 
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
 A logical argument indicating whether the iteration number should be printed. The default value is \code{FALSE}.
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
The \code{acenlm} function implements the ACE algorithm to find designs for general classes of nonlinear regression models with identically and independently normally distributed errors meaning the user does not have to write their own utility function.

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

In both cases the utility function is not available in closed form due to the analytical intractability of either the posterior distribution (for SIG) or the posterior mean (for NSEL). Sampling-based Monte Carlo or importance sampling approximations will be employed. This was the original approach used by Overstall & Woods (2017). 

A normal approximation to the posterior can be taken leading to the approximation by some scalar function of the Fisher information matrix, \eqn{\mathcal{I} (\theta;d)}{I (\theta;d)}, which only depends on \eqn{\theta}{\theta} (Chaloner & Verdinelli, 1995). In the case of SIG, the approximate utility is given by
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

Similar to all coordinate exchange algorithms, ACE should be repeated from different initial designs. The function 
\code{pacenlm} will implement this where the initial designs are given by a list via the argument \code{start.d}. On the completion 
of the repetitions of ACE, \code{pacenlm} will approximate the expected utility for all final designs and return the design (the terminal design) with the 
largest approximate expected utility.

For more details on the ACE algorithm, see Overstall & Woods (2017).
}

\value{
The function will return an object of class \code{"ace"} (for \code{acenlm}) or \code{"pace"} (for \code{pacenlm})  which is a list with the following components:
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
This will be \code{FALSE}.
}
\item{nlm}{
If the object is a result of a direct call to \code{acenlm} then this is \code{TRUE}.
}
\item{criterion}{
If the object is a result of a direct call to \code{acenlm} then this is the argument \code{criterion}.
}
\item{method}{
If the object is a result of a direct call to \code{acenlm} then this is the argument \code{method}.	
}
\item{prior}{
If the object is a result of a direct call to \code{aceglm} then this is the argument \code{prior}.
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
The terminal design (\code{pacenlm} only).
}
\item{eval}{
If \code{deterministic = "MC"}, a vector containing \code{n.assess} approximations to the expected utility for the terminal design (\code{pacenlm} only).

If \code{deterministic = "quadrature"}, a scalar giving the approximate expected utility for the terminal design (\code{pacenlm} only).
}
\item{final.d}{
A list of the same length as the argument \code{start.d}, where each element is the final design (i.e. \code{phase2.d}) for each repetition of the ACE algorithm (\code{pacenlm} only).
}
\item{besti}{
A scalar indicating which repetition of the ACE algorithm resulted in the terminal design (\code{pacenlm} only).
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
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\note{
These are wrapper functions for \code{\link{ace}} and \code{\link{pace}}.
}

\seealso{
\code{\link{ace}}, \code{\link{aceglm}}, \code{\link{pace}}, \code{\link{paceglm}}.
}
\examples{
## This example uses aceglm to find a Bayesian D-optimal design for a 
## compartmental model with 6 runs 1 factor. The priors are 
## those used by Overstall & Woods (2017). The design space for each 
## coordinate is [0, 24] hours.

set.seed(1)
## Set seed for reproducibility.

n<-6
## Specify the sample size (number of runs).

start.d<-matrix(24 * randomLHS(n = n, k = 1), nrow = n, ncol = 1,
dimnames = list(as.character(1:n), c("t")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.

low<-c(0.01884, 0.298, 21.8)
upp<-c(0.09884, 8.298, 21.8)
## Lower and upper limits of the support of the uniform prior distributions. Note that the prior
## for the third element is a point mass.

prior<-function(B){
out<-cbind(runif(n = B, min = low[1], max = upp[1]), runif(n = B, min = low[2],max = upp[2]),
runif(n = B, min = low[3], max = upp[3]))
colnames(out)<-c("a", "b", "c")
out}

## Create a function which specifies the prior. This function will return a 
## B by 3 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.

example1<-acenlm(formula = ~ c*(exp( - a * t) - exp( - b * t)), start.d = start.d, prior = prior, 
N1 = 1, N2 = 0, B = c(1000, 1000), lower = 0, upper = 24, method = "MC")
## Call the acenlm function which implements the ACE algorithm requesting 
## only one iteration of Phase I and zero iterations of Phase II. The Monte
## Carlo sample size for the comparison procedure (B[1]) is set to 1000.

example1
## Print out a short summary.

#Non Linear Model 
#Criterion = Bayesian D-optimality 
#Formula: ~c * (exp(-a * t) - exp(-b * t))
#Method:  MC 
#
#B:  1000 1000 
#
#Number of runs = 6
#
#Number of factors = 1
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 0
#
#Computer time = 00:00:00

example1$phase2.d
## Look at the final design.

#           t
#1 19.7787011
#2  2.6431912
#3  0.2356938
#4  8.2471451
#5  1.4742319
#6 12.7062270

prior2 <- list(support = cbind(rbind(low, upp)))
colnames(prior2$support) <- c("a", "b", "c")
example2 <- acenlm(formula = ~ c * (exp( - a * t) - exp( - b *t)), start.d = start.d, 
prior = prior2, lower = 0, upper = 24, N1 = 1, N2 = 0 )
## Call the acenlm function with the default method of "quadrature"

example2$phase2.d
## Final design

#           t
#1  0.5167335
#2  2.3194434
#3  1.5365409
#4  8.2471451
#5 21.9402670
#6 12.7062270

utility <- utilitynlm(formula = ~c * (exp( - a * t) - exp( - b *t)), prior = prior, 
                            desvars = "t", method = "MC" )$utility
## create a utility function to compare designs

mean(utility(example1$phase1.d, B = 20000))
#[1] 12.13773
mean(utility(example2$phase1.d, B = 20000))
#[1] 11.19745
## Compare the two designs using the Monte Carlo approximation
}
