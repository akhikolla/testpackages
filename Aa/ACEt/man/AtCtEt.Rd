\name{AtCtEt}
\alias{AtCtEt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting the ACE(t) model  
}
\description{
The ACE(t) model with the A, C and E variance components as functions with respect to age modelled by B-splines.
}
\usage{
AtCtEt(data_m, data_d, mod = c('d','d','d'), knot_a = 5, knot_c = 5, knot_e = 5, 
loc = c('e','e','e'), boot=FALSE, num_b = 100, init = rep(0,3), robust = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_m}{
An \eqn{N_m} x 3 data matrix for MZ twins. \eqn{N_m} is the number of MZ twin pairs. The first two columns are centered trait values (i.e. the mean should be zero) and the third column is age (or other covariates). 
}
  \item{data_d}{
An \eqn{N_d} x 3 data matrix for DZ twins. \eqn{N_d} is the number of DZ twin pairs. The first two columns are centered trait values (i.e. the mean should be zero) and the third column is age (or other covariates). 
}
  \item{mod}{
  A character vector of length 3. Each element specifies the function for the A, C or E component respectively. The A and C components can be 'd'(dynamic), 'c'(constant) or 'n'(NA). The E component can only be 'd' or 'c'. Thus, \eqn{model=c('c','c','c')} is corresponding to the classical ACE model.
  }
  \item{knot_a}{
The number of interior knots of the B-spline for the A component, which must be no less than 3. The default value is 5. 
}

  \item{knot_c}{
The number of interior knots of the B-spline for the C component, which must be no less than 3. The default value is 5.
}
  \item{knot_e}{
The number of interior knots of the B-spline for the E component, which must be no less than 3. The default value is 5.
}

  \item{loc}{
	A 1x3 character vector indicating how to place knots for each component: evenly ("e") or quantile-based ("q"). The default value is "e".
}

  \item{boot}{
A logical indicator of whether to use the bootstrap method to calculate the confidence interval. The default is FALSE. 
}

  \item{num_b}{
The number of replicates when the bootstrap method is used (i.e. \eqn{boot=TRUE}). The default value is 100.
}

 \item{init}{
A 3x1 vector of the initial values for the optimization. The default values are 1.
}

 \item{robust}{
An integer indicating the number of different initial values that the function will randomly generate and try in the optimization. The default value is 0.
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{n_beta_a }{The number of spline coefficients for the A component.}
  \item{n_beta_c }{The number of spline coefficients for the C component.}
  \item{n_beta_e }{The number of spline coefficients for the E component.}
  \item{beta_a }{The estimated spline coefficients (if the model parameter is 'd') or variance (if the model parameter is 'c') of the A component.}
  \item{beta_c }{The estimated spline coefficients (if the model parameter is 'd') or variance (if the model parameter is 'c') of the C component.}
  \item{beta_e }{The estimated spline coefficients (if the model parameter is 'd') or variance (if the model parameter is 'c') of the E component.}
  \item{hessian_ap }{The approximated Hessian matrix from the quasi-Newton algorithm.}
  \item{hessian }{The Hessian matrix calculated analytically.}
  \item{con }{An indicator of convergence of the optimization algorithm. An integer code 0 indicates successful completion. See 'optim' for more details.}
  \item{lik }{The minus log-likelihood.}
  \item{knots_a }{A vector of the knot positions for the A component.}
  \item{knots_c }{A vector of the knot positions for the C component.}
  \item{knots_e }{A vector of the knot positions for the E component.}
  \item{boot }{A list containing pointwise CIs estimated from the bootstrap method when \eqn{boot=TRUE}}
%% ...
}

\details{
If the variance is close to the boundary (0), it is better to use the bootstrap method to get the CIs.  
The optimization algorithm may sometimes end up with a local minimum. It is recommended to try different random initial values by setting 'robust'.
}


\references{
He, L., Sillanpää, M.J., Silventoinen, K., Kaprio, J. and Pitkäniemi, J., 2016. Estimating Modifying Effect of Age on Genetic and Environmental Variance Components in Twin Models. Genetics, 202(4), pp.1313-1328.

He, L., Pitkäniemi, J., Silventoinen, K. and Sillanpää, M.J., 2017. ACEt: An R package for estimating dynamic heritability and comparing twin models. Behavior Genetics, 47(6), pp.620-641.
}
\author{
Liang He
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(data_ace)

result <- AtCtEt(data_ace$mz, data_ace$dz, mod=c('d','d','c'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

