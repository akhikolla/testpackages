\name{size}

\alias{size}

\title{Computing the Size of Autocorrelation Robust Tests}
\encoding{latin1}
\description{
This function provides an implementation of Algorithm 2 in 
\enc{Pötscher}{Poetscher} and 
Preinerstorfer (2016), which we here abbreviate as (A2). The user is referred to
this article for definitions, a detailed description of the problem solved
by (A2), and for a detailed description of the algorithm itself.
}

\usage{
size(C, ar.order.max, bandwidth, ker, R, X, N0, N1, N2, Mp, M1, 
               M2, Eicker = FALSE, 
               opt.method.1 = "Nelder-Mead",
               opt.method.2 = "Nelder-Mead", 
               control.1 = list("reltol" = N1^(-.5), "maxit" = dim(X)[1]*20),
               control.2 = list("reltol" = N2^(-.5), "maxit" = dim(X)[1]*30),   
               cores = 1, margin = rep(1, length = ar.order.max))
}

\arguments{
\item{C}{Critical value. A positive real number (for negative critical values
the size of the test equals \eqn{1}).}
                 
\item{ar.order.max}{Maximal order of the stationary autoregressive 
error process. A nonnegative integer. 
If \code{ar.order.max} is greater than sample size minus one, then 
\code{ar.order.max} is 
automatically replaced by sample size minus one. If \code{ar.order.max} equals
sample size minus one, then the size approximated by (A2) is with respect to 
the set of all stationary Gaussian error processes. If \code{ar.order.max} is 
set equal to 0, \code{size} does not return the outcome of (A2), 
because no optimization is needed in that case. Instead, 
it draws a pseudorandom sample of size \code{N2} from the distribution of the 
test statistic under the null hypothesis and returns the fraction of
observations that give rise to a test statistic not smaller than the critical 
value \code{C}.} 
         
\item{bandwidth}{Bandwidth parameter used in the construction of the test 
statistic. A positive real number.}

\item{ker}{Kernel function used in the construction of the test statistic. 
\code{ker} can take one of the
values "Bartlett", "Parzen", or "Quadratic Spectral". The 
\code{\link{kweights}}
function is used to generate 
the weights. For 
the test statistic \eqn{T_{w}} this implies the weights used via 
\eqn{w(j, n) = ker(j/bandwidth)}. For  the test statistic 
\eqn{T_{E, \mathsf{W}}} this
implies the weights matrix \eqn{\mathsf{W}} via \eqn{\mathsf{W}_{ij} = 
ker((i-j)/bandwidth)}, cf. 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016) for definitions of
\eqn{T_{w}} and
\eqn{T_{E, \mathsf{W}}}.}

\item{R}{The restriction matrix. \code{size} computes the size of a test 
for the hypothesis \eqn{R \beta = r}. \code{R} needs to be of
full row rank, and needs to have the same number of columns as \code{X}.}
                 
\item{X}{The design matrix. \code{X} needs to be of full column rank. 
The number of columns of \code{X} must be smaller than the number of rows of 
\code{X}.}

\item{N0}{A positive integer. Corresponds to \eqn{N_0} in the description of 
(A2) in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016).}

\item{N1}{A positive integer. Corresponds to \eqn{N_1} in the description of 
(A2) in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016). N1 should be greater than 
N0.}     

\item{N2}{A positive integer. Corresponds to \eqn{N_2} in the description of 
(A2) in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016). N2 should be greater than
N1.} 

\item{Mp}{A positive integer. \code{Mp} determines \eqn{M_0} in (A2), i.e., 
the number of starting 
values 
chosen in 
Stage 0 of (A2). The way initial values are generated depends on 
\code{ar.order.max}: 
\enumerate{
\item If \code{ar.order.max} is \eqn{0} the choice of \code{Mp} is irrelevant, 
see the
description of \code{ar.order.max} above. 
\item If \code{ar.order.max} is either \eqn{1} or \eqn{2}, then the initial 
values 
are an i.i.d. pseudorandom sample of size \code{Mp} drawn from a distribution on
the 
set of partial autocorrelation coefficients that induces a uniform distribution
on the set of stationary autoregressive coefficients of order 
\code{ar.order.max} (cf. Jones (1987)). 
\item If \code{ar.order.max} is greater than \eqn{2}, then starting values are 
generated as follows: Let \eqn{A} denote the set of integer
multiples of \eqn{5} that are
greater than \eqn{2} and smaller than \code{ar.order.max}.
For every number \deqn{l \in A \cup \{2, ar.order.max\}} the algorithm in Jones 
(1987) is used to generate a pseudorandom sample of size \code{Mp} on the set of
partial autocorrelation coefficients from a distribution that induces a uniform
distribution on the set of stationary AR(\eqn{l}) coefficients.
These samples are then
used as initial values, cf. also the discussion in \enc{Pötscher}{Poetscher} and 
Preinerstorfer (2016).}
}

\item{M1}{A positive integer. Corresponds to \eqn{M_1} in the description of 
(A2) in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016). M1 must not exceed Mp.}

\item{M2}{A positive integer. Corresponds to \eqn{M_2} in the description of 
(A2) in 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016). M2 must not exceed M1.}

\item{Eicker}{Determines the test for which the size is computed. If 
\code{Eicker = TRUE}, then 
the size of the test that rejects if \eqn{T_{E, \mathsf{W}}} is not smaller than
\code{C} (with 
\eqn{\mathsf{W}_{ij}=ker((i-j)/bandwidth)}) is computed. If 
\code{Eicker = FALSE}, then 
the size of the test that rejects if \eqn{T_{w}} is not smaller than
\code{C} (with \eqn{w(j, n) = 
ker(j/bandwidth)}) is computed (cf. \enc{Pötscher}{Poetscher} and
Preinerstorfer 
(2016) for a precise definition of these test statistics). Default is 
\code{Eicker = FALSE}.} 
    
\item{opt.method.1}{The optimization method chosen in Stage 1
of (A2). Any optimization routine implemented in \code{\link{optim}} apart from
"Brent" can be supplied. The default is "Nelder-Mead".} 

\item{opt.method.2}{The optimization method chosen in Stage 2
of (A2). Any optimization routine implemented in \code{\link{optim}} apart from
"Brent" can be supplied. The default is "Nelder-Mead".} 
 
\item{control.1}{Control parameters passed to the \code{optim} function in 
Stage 1 of (A2). Default is \code{control.1 = list("reltol" = N1^(-.5), 
"maxit" = dim(X)[1]*20)}.}

\item{control.2}{Control parameters passed to the \code{optim} function in 
Stage 2 of (A2). Default is \code{control.2 = list("reltol" = N2^(-.5), 
"maxit" = dim(X)[1]*30)}.}
 
\item{cores}{The number of CPU cores used in the (parallelized)
computation of the Monte-Carlo
approximations in (A2). Default is 1. Parallelized computation is enabled only
if the compiler used to build \pkg{acrt} supports OpenMP.}
 
\item{margin}{The restrictions imposed on the partial autocorrelation 
coefficients. \code{margin} is an \code{ar.order.max}-dimensional vector of 
real numbers in 
\eqn{(0, 1]}. Default is \code{margin = rep(1, length = ar.order.max)}, 
which corresponds to no restriction.}

}

\details{For details see the relevant sections in \enc{Pötscher}{Poetscher} 
and Preinerstorfer 
(2016), in particular the description of Algorithm 2 in the Appendix.
}

\value{The output of \code{size} depends on \code{ar.order.max}. If
\code{ar.order.max} is zero, the function \code{size}
returns a list consisting of:

\item{size}{The size obtained by drawing a pseudorandom sample of size \code{N2}
from the distribution of the test statistic under the null hypothesis and 
by computing the fraction of
observations giving rise to a test statistic not smaller than the critical value
\code{C}.
}

If \code{ar.order.max} is greater than zero, the function \code{size}
returns a list consisting of:

\item{starting.parameters}{The rows of this matrix are the
initial values (partial autocorrelation 
coefficient vectors) that were used in Stage 1 of (A2), and which were 
chosen from the pool of randomly generated initial values in Stage 0. The rows
correspond to \eqn{\rho_{1:M_0}, ..., \rho_{M_1:M_0}}, respectively, in the 
description of (A2).}

\item{starting.rejection.probs}{Monte Carlo approximations of the null-rejection
probabilities corresponding to the initial values used in Stage 1 of (A2).
The coordinates
of this vector correspond to 
\eqn{\tilde{p}_{1:M_0}, ..., \tilde{p}_{M_1:M_0}} 
in the description of (A2).}

\item{first.stage.parameters}{The rows of this matrix are the 
parameters (partial autocorrelation 
coefficients) that were obtained in Stage 1 of (A2). The rows
correspond to \eqn{\rho^*_{1}, ..., \rho^*_{M_1}}, respectively, in the 
description of (A2).}

\item{first.stage.rejection.probs}{Monte Carlo approximations of the 
null-rejection
probabilities corresponding to the \code{first.stage.parameters}, i.e.,
\eqn{\bar{p}_{j, \rho_j^*}} 
for \eqn{j = 1, ..., M_1} in the description of (A2).} 

\item{second.stage.parameters}{The rows of this matrix are the 
parameters (partial autocorrelation 
coefficients) that were obtained in Stage 2 of (A2). 
The rows correspond to \eqn{\rho^{**}_{1}, ..., \rho^{**}_{M_2}}, respectively,
in the description of (A2).}

\item{second.stage.rejection.probs}{Monte Carlo approximations of the null
rejection probabilities
corresponding to the \code{second.stage.parameters}, i.e.,
\eqn{\bar{\bar{p}}_{j, \rho_j^{**}}} for \eqn{
j = 1, ..., M_2} in the description of (A2).}

\item{convergence}{Convergence codes returned from \code{\link{optim}}
in Stage 2 of (A2) for each initial value.}

\item{size}{The size computed by (A2), i.e., the 
maximum of \eqn{\bar{\bar{p}}_{j, \rho_j^{**}}} for 
\eqn{j = 1, ..., M_2}.}
}

\examples{
\dontrun{
n <- 100
C <- 2.260568^2
ar.order.max <- n-1
bandwidth <- n/10
ker <- "Bartlett"
R <- matrix(c(0, 1), nrow = 1, ncol = 2)
X <- cbind(rep(1, length = n), rnorm(n))
N0 <- 1000
N1 <- 10000
N2 <- 50000
Mp <- 5000
M1 <- 10
M2 <- 2

size(C, ar.order.max, bandwidth, ker, R, X, N0, N1, N2, Mp, M1, M2)
}}

\references{
Jones, M. C. (1987). Randomly choosing parameters from the stationarity and 
invertibility region of autoregressive-moving average models. \emph{Applied 
Statistics}, \strong{36} 134-138

\enc{Pötscher}{Poetscher}, B.M. and Preinerstorfer, D. (2016). Controlling
the size of autocorrelation robust tests.
\url{https://arxiv.org/abs/1612.06127/}
}

\seealso{
\code{\link{optim}}, \code{\link{kweights}}.
}                  