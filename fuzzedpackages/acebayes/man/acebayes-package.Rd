\name{acebayes-package}
\alias{acebayes-package}
\alias{acebayes}
\docType{package}
\title{
Optimal Bayesian Experimental Design using the Approximate Coordinate Exchange (ACE) Algorithm
}
\description{
Finding an optimal Bayesian experimental design (Chaloner & Verdinelli, 1995) involves maximising 
an objective function given by the expectation of some appropriately chosen utility 
function with respect to the joint distribution of unknown quantities (including responses). This objective 
function is usually not available in closed form and the design space can be continuous and of high 
dimensionality.

The \code{acebayes} package uses Approximate Coordinate Exchange (ACE; Overstall 
& Woods, 2017) to maximise an approximation to the expectation of the utility function. 
In Phase I of the algorithm, a continuous version of the coordinate exchange 
algorithm (Meyer & Nachtsheim, 1995) is used to maximise an approximation to 
expected utility. The approximation is given by the predictive mean of a Gaussian 
process (GP) emulator constructing using a 'small' number of approximate
evaluations of the expected utility function. In Phase II a point exchange 
algorithm is used to consolidate clusters of design points into repeated design 
points.

}
\details{
\tabular{ll}{
Package: \tab acebayes\cr
Version: \tab 1.10\cr
Date: \tab 2020-10-04\cr
Date: \tab 2017-02-09\cr
License: \tab GPL-2\cr
}
The most important functions are as follows.

\enumerate{
\item
\code{\link{ace}}
\item
\code{\link{pace}}
\item
\code{\link{aceglm}}
\item
\code{\link{paceglm}}
\item
\code{\link{acenlm}}
\item
\code{\link{pacenlm}}
}

The function \code{\link{ace}} implements both phases of the ACE algorithm. It has two mandatory arguments: \code{utility} (a function specifying the chosen utility function incorporating the joint distribution of unknown quantities) and \code{start.d} (the initial design). The function will return the final design from the algorithm, along with information to assess convergence. The function \code{\link{pace}} implements repetitions of the ACE algorithm from different starting designs (as specified by the \code{start.d} argument).

The computational time of \code{\link{ace}} (and \code{\link{pace}}) is highly dependent on the computational time required to evaluate the user-supplied function \code{utility}. Therefore it is recommended that users take advantage of \code{R} packages such as \code{Rcpp} (Eddelbuettel & Francois, 2011), \code{RcppArmadillo} (Eddelbuettel & Sanderson, 2014), or \code{RcppEigen} (Bates & Eddelbuettel, 2013), that provide convenient interfaces to compiled programming languages.

The functions \code{\link{aceglm}} and \code{\link{acenlm}} are user-friendly wrapper functions for \code{\link{ace}} which use the ACE algorithm to find Bayesian optimal experimental designs for generalised linear models and non-linear models, respectively. As special cases, both of these functions can find pseudo-Bayesian optimal designs. The functions \code{\link{paceglm}} and \code{\link{pacenlm}} implement repetitions of the ACE algorithm from different starting designs (as specified by the \code{start.d} argument) for generalised linear models and non-linear models, respectively.

For more details on the underpinning methodology, see Overstall & Woods (2017), and for more information on the \code{acebayes} package, see Overstall et al (2020).
}
\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides

Maintainer: Antony. M.Overstall \email{A.M.Overstall@soton.ac.uk}
}
\references{
Bates, D. & Eddelbuettel, D. (2013). Fast and Elegant Numerical Linear Algebra Using the 
RcppEigen Package. \emph{Journal of Statistical Software}, \bold{52(5)}, 1-24. 
\url{https://www.jstatsoft.org/v52/i05/}

Chaloner, K. & Verdinelli, I. (1995). Bayesian Experimental Design: A Review. \emph{Statistical
Science}, \bold{10}, 273-304.

Eddelbuettel, D. & Francois, R. (2011). Rcpp: Seamless R and C++ Integration. 
\emph{Journal of Statistical Software}, \bold{40(8)}, 1-18. \url{https://www.jstatsoft.org/v40/i08/}

Eddelbuettel, D. & Sanderson, C. (2014). RcppArmadillo: Accelerating R with high-performance 
C++ linear algebra. \emph{Computational Statistics and Data Analysis}, \bold{71}, 1054-1063.

Meyer, R. & Nachtsheim, C. (1995). The Coordinate Exchange Algorithm for Constructing Exact
Optimal Experimental Designs. \emph{Technometrics}, \bold{37}, 60-69.

Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.

Overstall, A.M., Woods, D.C. & Adamou, M. (2020). acebayes: An R Package for Bayesian Optimal Design of Experiments via Approximate Coordinate Exchange. \emph{Journal of Statistical Software}, \bold{95} (13), 1-33 \url{https://www.jstatsoft.org/v095/i13/}
}

\keyword{package}

\examples{
## This example uses aceglm to find a pseudo-Bayesian D-optimal design for a 
## first-order logistic regression model with 6 runs 4 factors (i.e. 5 parameters).
## The priors are those used by Overstall & Woods (2017), i.e. a uniform prior 
## distribution is assumed for each parameter. The design space for each coordinate 
## is [-1, 1].

set.seed(1)
## Set seed for reproducibility.

n<-6
## Specify the sample size (number of runs).

start.d<-matrix(2 * randomLHS(n = n,k = 4) - 1, nrow = n, ncol = 4,
dimnames = list(as.character(1:n), c("x1", "x2", "x3", "x4")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.

low<-c(-3, 4, 5, -6, -2.5)
upp<-c(3, 10, 11, 0, 3.5)
## Lower and upper limits of the uniform prior distributions.

prior<-function(B){
t(t(6*matrix(runif(n = 5*B), ncol = 5)) + low)}
## Create a function which specifies the prior. This function will return a 
## B by 5 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.

example<-aceglm(formula = ~ x1 + x2 + x3 + x4,  start.d = start.d, family = binomial, 
prior = prior , criterion = "D", method= "MC", B = c(1000,1000), N1 = 1, N2 = 0, 
upper = 1)
## Call the aceglm function which implements the ACE algorithm requesting 
## only one iteration of Phase I and zero iterations of Phase II (chosen for 
## illustrative purposes). The Monte Carlo sample size for the comparison 
## procedure (B[1]) is set to 1000 (chosen again for illustrative purposes).

example
## Print out a short summary. 

#Generalised Linear Model
#Criterion = Bayesian D-optimality 
#
#Number of runs = 6
#
#Number of factors = 4
#
#Number of Phase I iterations = 1
#
#Number of Phase II iterations = 0
#
#Computer time = 00:00:02

## The final design found is:

example$phase2.d

#          x1          x2          x3         x4
#1 -0.4735783  0.12870470 -0.75064318  1.0000000
#2 -0.7546841  0.78864527  0.58689270  0.2946728
#3 -0.7463834  0.33548985 -0.93497463 -0.9573198
#4  0.4446617 -0.29735212  0.74040030  0.2182800
#5  0.8459424 -0.41734194 -0.07235575 -0.4823212
#6  0.6731941  0.05742842  1.00000000 -0.1742566
}
