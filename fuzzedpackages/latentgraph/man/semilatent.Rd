\name{semilatent}
\alias{semilatent}

\title{
Estimate Graphical Models with Latent Variables and Replicates
}
\description{
Estimate graphical models with latent variables and replicates using the method in Tan et al. (2016).
}
\usage{
semilatent(data, n, R, p, lambda, distribution = "Gaussian", rule = "AND")
}

\arguments{
  \item{data}{data set. Can be a matrix, list, array, or data frame. If the data set is a matrix, it should have \eqn{nR} rows and \eqn{p} columns. This matrix is formed by stacking \eqn{n} matrices, and each matrix has \eqn{R} rows and \eqn{p} columns. If the data set is a data frame, the dimention and structure are the same as the matrix. If the data set is an array, its dimention is (R, p, n). If the data set is a list, it should have \eqn{n} elements and each element is a matrix with \eqn{R} rows and \eqn{p} columns.}
  \item{n}{the number of observations.}
  \item{R}{the number of replicates for each observation.}
  \item{p}{the number of observed variables.}
  \item{lambda}{tuning parameter that encourages estimated graph to be sparse.}
  \item{distribution}{For a data set with Gaussian distribution, use "Gaussian"; For a data set with Ising distribution, use "Ising". Default is "Gaussian".}
  \item{rule}{rules to combine matrices that encode the conditional dependence relationships between sets of two observed variables. Options are "AND" and "OR". Default is "AND".}
}

\details{
The semilatent method has two assumptions. The first one states that the latent variables are constant across replicates.
Assumption 2 states that given the latent variables, the replicates are mutually independent.
With these two assumptions, the method seeks to solve the following problem for \eqn{1 \le j \le p}.
\deqn{
\min_{\beta_{j,O / j}} \{l_j (\beta_{j,O / j}) + \lambda\|\beta_{j,O / j}\|_1 \},
}
where \eqn{l_j (\beta_{j,O / j}) } is a nuisance-free loss function, \eqn{\beta_{j,O / j}} is a parameter that represents the conditional dependence relationships between \eqn{j}th observed variable and the other observed variables, and \eqn{\lambda} is a tuning parameter.
This method aims at modeling semiparametric exponential family graphical model with latent variables and replicates.
}

\value{
  \item{omega}{a matrix that encodes the conditional dependence relationships between sets of two observed variables}
  \item{theta}{the adjacency matrix with 0 and 1 encoding conditional independence and dependence between sets of two observed variables, respectively}
  \item{penalty}{the penalty value}
}

\references{
Tan, K. M., Ning, Y., Witten, D. M. & Liu, H. (2016), `Replicates in high dimensions, with applications to latent variable graphical models', Biometrika 103(4), 761--777.
}


\examples{
#semilatent Gaussian with "AND" rule
n <- 50
R <- 20
p <- 30
seed <- 1
l <- 2
s <- 2
data <- generate_Gaussian(n, R, p, l, s, sparsityA = 0.95, sparsityobserved = 0.9,
sparsitylatent = 0.2, lwb = 0.3, upb = 0.3, seed)$X

result <- semilatent(data, n, R, p, lambda = 0.1,distribution = "Gaussian", 
rule = "AND")
}

