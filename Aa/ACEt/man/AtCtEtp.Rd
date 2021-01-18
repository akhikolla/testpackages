\name{AtCtEtp}
\alias{AtCtEtp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting the ACE(t)-p model 
}
\description{
The ACE(t)-p model with the A, C and E variance components as functions with respect to age modelled by P-splines.
}
\usage{
AtCtEtp(data_m, data_d, knot_a = 8, knot_c = 8, knot_e = 8, eps = 0.1, 
mod=c('d','d','d'), robust = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_m}{
An \eqn{N_m} x 3 data matrix for MZ twins. \eqn{N_m} is the number of MZ twin pairs. The first two columns are centered trait values (i.e. the mean should be zero) and the third column is age (or other covariates). 

}
  \item{data_d}{
An \eqn{N_d} x 3 data matrix for DZ twins. \eqn{N_d} is the number of DZ twin pairs. The first two columns are centered trait values (i.e. the mean should be zero) and the third column is age (or other covariates). 
}
  \item{knot_a}{
The number of interior knots of the B-spline for the A component. The default value is 8.
}
  \item{knot_c}{
The number of interior knots of the B-spline for the C component. The default value is 8.
}

  \item{knot_e}{
The number of interior knots of the B-spline for the E component. The default value is 8.
}
  \item{eps}{
Tolerance for convergence of the EM algorithm iterations. The default value is 0.1.
}

  \item{mod}{
A character vector of length 3. Each element specifies the function for the A, C or E component respectively. The components can be 'd'(dynamic), 'c'(constant) or 'l'(linear). The default is c('d','d','d').
}
  \item{robust}{
  An integer indicating the number of different initial values that the function will randomly generate and try in the optimization. The default value is 0.
  }

}

\value{
  
  \item{var_b_a }{The estimated variance for the penalized coefficient for the A components.}
  \item{var_b_c }{The estimated variance for the penalized coefficient for the C components.}
  \item{var_b_e }{The estimated variance for the penalized coefficient for the E components.}
  \item{beta_a }{The estimated spline coefficients of the A component. See 'details' for more information.}
  \item{beta_c }{The estimated spline coefficients of the C component. See 'details' for more information.}
  \item{beta_e }{The estimated spline coefficients of the E component. See 'details' for more information.}
  \item{con }{An indicator of convergence of the optimization algorithm. An integer code 0 indicates successful completion. See 'optim' for more details.}
  \item{lik }{The minus log marginal likelihood.}
  \item{knot_a }{A vector of the knot positions for the A component.}
  \item{knot_c }{A vector of the knot positions for the C component.}
  \item{knot_e }{A vector of the knot positions for the E component.}
%% ...
}

\details{
When the 'mod' argument for a component is 'd'(dynamic), the corresponding 'beta' is the spline coefficients.   
When the 'mod' argument for a component is 'l'(linear), the corresponding 'beta' is a vector of two values, the exponential of which (exp(beta)) are the variances at the minimum and maximum age (or other covariates) provided in the data.     
When the 'mod' argument for a component is 'c'(constant), the corresponding 'beta' has only one value and exp(beta) is the variance.
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
# data(data_ace)

# result <- AtCtEtp(data_ace$mz, data_ace$dz, knot_e = 7, knot_c = 5, mod=c('d','d','d'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

