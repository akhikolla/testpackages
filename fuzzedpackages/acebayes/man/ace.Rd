\name{ace}
\alias{ace}
\alias{pace}
\alias{acephase1}
\alias{acephase2}

\title{
Approximate Coordinate Exchange (ACE) Algorithm
}
\encoding{UTF-8}
\description{
These functions implement the approximate coordinate exchange (ACE) algorithm (Overstall & Woods, 2017) for finding optimal Bayesian experimental designs by maximising an approximation to an intractable expected utility function.
}
\usage{
ace(utility, start.d, B, Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
limits = NULL, progress = FALSE, binary = FALSE, deterministic = FALSE)

acephase1(utility, start.d, B, Q = 20, N1 = 20, lower, upper, limits = NULL, 
progress = FALSE, binary = FALSE, deterministic = FALSE)

acephase2(utility, start.d, B, N2 = 100, progress = FALSE, binary = FALSE, 
deterministic = FALSE)

pace(utility, start.d, B, Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, 
limits = NULL, binary = FALSE, deterministic = FALSE, mc.cores = 1, 
n.assess = 20)
}

\arguments{
  \item{utility}{
A function with two arguments: \code{d} and \code{B}. 

For a Monte Carlo approximation (\code{deterministic = FALSE}), it should return a \emph{\bold{vector}} of length \code{B} where each element gives the value of the utility function for design \code{d}, for a value generated from the joint distribution of all unknown quantities. The mean of the elements of this vector provides a Monte Carlo approximation to the expected utility for design \code{d}.

For a deterministic approximation (\code{deterministic = TRUE}), it should return a \emph{\bold{scalar}} giving the approximate value of the expected utility for design \code{d}. In this latter case, the argument \code{B} can be a list containing tuning parameters for the deterministic approximation. If \code{B} is not required, the utility function must still accept the argument.
}
  \item{start.d}{
For \code{ace}, \code{acephase1} and \code{acephase2}, an \code{n} by \code{k} matrix specifying the initial design for the ACE algorithm.

For \code{pace}, a list with each element being an \code{n} by \code{k} matrix specifying the initial design for each repetition of the ACE algorithm.
}
  \item{B}{
An argument for controlling the approximation to the expected utility.

For a Monte Carlo approximation (\code{deterministic = FALSE}), a vector of length two specifying the size of the Monte Carlo samples, generated from the joint distribution of unknown quantities. The first sample size, \code{B[1]}, gives the sample size to use in the comparison procedures, and the second sample size, \code{B[2]}, gives the sample size to use for the evaluations of Monte Carlo integration that are used to fit the Gaussian process emulator. If missing when \code{deterministic = FALSE}, the default value is \code{c(20000,1000)}. 

For a deterministic approximation (\code{deterministic = TRUE}), then \code{B} may be a list of length two containing any necessary tuning parameters for the expected utility calculations for the comparison and emulation steps.
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
An argument specifying the bounds on the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the lower limits of all coordinates of the design space. The default value is -1.
}
  \item{upper}{
An argument specifying the bounds on the design space. This argument can either be a scalar or a matrix of the same dimension as the argument \code{start.d} which specifies the upper limits of all coordinates of the design space. The default value is 1.
}
  \item{limits}{
An argument specifying the grid over which to maximise the Gaussian process emulator for the expected utility function. It should be a function with three arguments: \code{i}, \code{j} and \code{d} which generates a one-dimensional grid for the \code{ij}th coordinate of the design when the current design is \code{d}. The default value is NULL which generates values uniformly on the interval \code{(lower[i,j],upper[i,j])} or \code{(lower,upper)} depending on whether the arguments \code{lower} and \code{upper} are matrices or scalars, respectively.
}
  \item{progress}{
A logical argument indicating whether the iteration number and other information detailing the progress of the algorithm should be printed. The default value is \code{FALSE}.
}
  \item{binary}{
A logical argument indicating whether the utility function has binary or continuous output. In some cases, the utility function is an indicator function of some event giving binary output. The expected utility function will then be the expected posterior probability of the event. Utility functions such as Shannon information gain and negative squared error loss give continuous output. The type of output guides the choice of comparison procedure used in the ACE algorithm. The default value is \code{FALSE}, indicating the utility function has continuous output.
}
	\item{deterministic}{
A logical argument indicating if a Monte Carlo (\code{FALSE}, default) or deterministic (\code{TRUE}) approximation to the expected utility is being used. 		
}
	\item{mc.cores}{
The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one (the default), and parallelisation requires at least two cores. See \code{\link[parallel]{mclapply}} for more information and warnings for \code{mc.cores} > 1. 	
}
	\item{n.assess}{
If \code{deterministic = FALSE}, the approximate expected utility for the design from each repetition of the ACE algorithm will be calculated \code{n.assess} times. The terminal design returned will be the design with the largest mean approximate expected utility calculated over the \code{n.assess} approximations.}
	}

\details{
Finding an optimal Bayesian experimental design (Chaloner & Verdinelli, 1995) involves maximising 
an objective function given by the expectation of some appropriately chosen utility 
function with respect to the joint distribution of unknown quantities (including responses). This objective 
function is usually not available in closed form and the design space can be continuous and of high 
dimensionality.

Overstall & Woods (2017) proposed the approximate coordinate exchange (ACE) algorithm 
to approximately maximise the expectation of the utility function. ACE consists of
two phases.

Phase I uses a continuous version of the coordinate exchange algorithm (Meyer & 
Nachtsheim, 1995) to maximise an approximation to the expected utility. Very briefly, 
the approximate expected utility is sequentially maximised over each one-dimensional element
of the design space. The approximate expected utility is given by the predictive mean of a 
Gaussian process (GP) regression model (also known as an emulator or surrogate) fitted 
to a 'small' number (argument \code{Q}) of evaluations of either a Monte Carlo (MC) or deterministic (e.g. quadrature) approximation to the expected utility (the MC sample size or arguments for the deterministic approximation are given by \code{B}). A GP 
emulator is a statistical model and, similar to all statistical models, can be an 
inadequate representation of the underlying process (i.e. the expected utility). 
Instead of automatically accepting the new design given by the value that maximises 
the GP emulator, for MC approximations a Bayesian hypothesis test, independent of the GP emulator, is 
performed to assess whether the expected utility of the new design is larger than the 
current design. For deterministic approximations, the approximate expected utility is calculated for the new design, and compared to that for the current design.

Phase I tends to produce clusters of design points. This is where, for example, two
design points are separated by small Euclidean distance. Phase II allows these points 
to be consolidated into a single repeated design point by using a point exchange 
algorithm (e.g Gotwalt et al., 2009) with a candidate set given by the final design 
from Phase I. In the same way as Phase I, comparisons of the expected loss between 
two designs is made on the basis of either a Bayesian hypothesis test or a direct comparison of deterministic approximations.

The original Bayesian hypothesis test proposed by Overstall & Woods (2017) is 
appropriate for utility functions with continuous output. Overstall et al. (2017) 
extended the idea to utility functions with binary output, e.g. the utility 
function is an indicator function for some event. The type of test can be specified by
the argument \code{binary}.

Similar to all coordinate exchange algorithms, ACE should be repeated from different initial designs. The function 
\code{pace} will implement this where the initial designs are given by a list via the argument \code{start.d}. On the completion 
of the repetitions of ACE, \code{pace} will approximate the expected utility for all final designs and return the design (the terminal design) with the 
largest approximate expected utility.
}

\value{
The function will return an object of class \code{"ace"} (for functions \code{ace}, \code{acephase1} and \code{acephase2}) or \code{"pace"} (for function \code{"pace"}) which is a list with the following components:

\item{utility}{
The argument \code{utility}.
}
\item{start.d}{
The argument \code{start.d}.
}
\item{phase1.d}{
The design found from Phase I of the ACE algorithm (only for \code{ace}, \code{acephase1} and \code{acephase2}).
}
\item{phase2.d}{
The design found from Phase II of the ACE algorithm (only for \code{ace}, \code{acephase1} and \code{acephase2})..
}
\item{phase1.trace}{
A vector containing the approximated expected utility of the current design at each stage of Phase I of the ACE algorithm. This can be used to assess convergence for MC approximations.
If \code{deterministic = FALSE}, this will be the mean of a call to \code{utility} with \code{d} being the current design and \code{B} being equal to the argument \code{B[1]}. If \code{deterministic = TRUE}, 
this will be a call to \code{utility} with \code{d} being the current design.

For \code{pace}, this will be \code{phase1.trace} for the terminal design.
}
\item{phase2.trace}{
A vector containing the approximated expected utility of the current design at each stage of Phase II of the ACE algorithm. This can be used to assess convergence for MC approximations.
If \code{deterministic = FALSE}, this will be the mean of a call to \code{utility} with \code{d} being the current design and \code{B} being equal to the argument \code{B[1]}. If \code{deterministic = TRUE}, 
this will be a call to \code{utility} with \code{d} being the current design.

For \code{pace}, this will be \code{phase2.trace} for the terminal design.
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
If the object is a result of a direct call to \code{ace} then this is \code{FALSE}.
}
\item{nlm}{
If the object is a result of a direct call to \code{ace} then this is \code{FALSE}.
}
\item{criterion}{
If the object is a result of a direct call to \code{ace} then this is \code{"NA"}.
}
\item{prior}{
If the object is a result of a direct call to \code{ace} then this is \code{"NA"}.
}
\item{time}{
Computational time (in seconds) to run the ACE algorithm.
}
\item{binary}{
The argument \code{binary}.
}
\item{deterministic}{
The argument \code{deterministic}.	
}
\item{d}{
The terminal design (\code{pace} only).
}
\item{eval}{
If \code{deterministic = FALSE}, a vector containing \code{n.assess} approximations to the expected utility for the terminal design (\code{pace} only).

If \code{deterministic = TRUE}, a scalar giving the approximate expected utility for the terminal design (\code{pace} only).
}
\item{final.d}{
A list of the same length as the argument \code{start.d}, where each element is the final design (i.e. \code{phase2.d}) for each repetition of the ACE algorithm (\code{pace} only).
}
\item{besti}{
A scalar indicating which repetition of the ACE algorithm resulted in the terminal design (\code{pace} only).
}
}
\references{
Chaloner, K. & Verdinelli, I. (1995). Bayesian Experimental Design: A Review. \emph{Statistical
Science}, \bold{10}, 273-304.

Gotwalt, C., Jones, B. & Steinberg, D. (2009). Fast Computation of Designs Robust to Parameter
Uncertainty for Nonlinear Settings. \emph{Technometrics}, \bold{51}, 88-95.

Meyer, R. & Nachtsheim, C. (1995). The Coordinate Exchange Algorithm for Constructing Exact
Optimal Experimental Designs. \emph{Technometrics}, \bold{37}, 60-69.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.

Overstall, A.M., McGree, J.M. & Drovandi, C.C. (2018). An approach for finding fully Bayesian optimal designs using normal-based approximations to loss functions. \emph{Statistics and Computing}, \bold{28}(2), 343-358.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\note{
For more details on the \code{R} implementation of the utility function used in the \bold{Examples} section, see \code{\link{utilcomp18bad}}.
}

\examples{
set.seed(1)
## Set seed for reproducibility.

## This example involves finding a pseudo-Bayesian D-optimal design for a 
## compartmental model with n = 18 runs. There are three parameters. 
## Two parameters have uniform priors and the third has a prior 
## point mass. For more details see Overstall & Woods (2017).

start.d<-optimumLHS(n = 18, k = 1)
## Create an initial design.

## Using a MC approximation
example1<-ace(utility = utilcomp18bad, start.d = start.d, N1 = 1, N2 = 2, B = c(100, 20))
## Implement the ACE algorithm with 1 Phase I iterations and 2 Phase II
## iterations. The Monte Carlo sample sizes are 100 (for comparison) and 20 for
## fitting the GP emulator.

example1
## Produce a short summary.

#User-defined model & utility 
#
#Number of runs = 18
#
#Number of factors = 1
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 2
#
#Computer time = 00:00:00

mean(utilcomp18bad(d = example1$phase2.d, B = 100))
## Calculate an approximation to the expected utility for the final design.
## Should get:

#[1] 9.254198

\dontrun{
plot(example1)
## Produces a trace plot of the current value of the expected utility. This
## can be used to assess convergence.
}

}