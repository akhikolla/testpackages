\name{corlatent}
\alias{corlatent}

\title{
Graphical Models with Latent Variables and Correlated Replicates
}
\description{
Estimate graphical models with latent variables and correlated replicates using the method in Jin et al. (2020).
}
\usage{
corlatent(data, accuracy, n, R, p, lambda1, lambda2, lambda3, distribution = "Gaussian",
rule = "AND")
}

\arguments{
  \item{data}{data set. Can be a matrix, list, array, or data frame. If the data set is a matrix, it should have \eqn{nR} rows and \eqn{p} columns. This matrix is formed by stacking \eqn{n} matrices, and each matrix has \eqn{R} rows and \eqn{p} columns. If the data set is a data frame, the dimention and structure are the same as the matrix. If the data set is an array, its dimention is (R, p, n). If the data set is a list, it should have \eqn{n} elements and each element is a matrix with \eqn{R} rows and \eqn{p} columns.}
  \item{accuracy}{the threshhold where algorithm stops. The algorithm stops when the difference between estimaters of the \eqn{(k-1)}th iteration and the \eqn{k}th iteration is smaller than the value of accuracy.}
  \item{n}{the number of observations.}
  \item{R}{the number of replicates for each observation.}
  \item{p}{the number of observed variables.}
  \item{lambda1}{tuning parameter that encourages estimated graph to be sparse.}
  \item{lambda2}{tuning parameter that models the effects of correlated replicates. Usually set to be equal to lambda1.}
  \item{lambda3}{tuning parameter that encourages the latent effect to be piecewise constants.}
  \item{distribution}{For a data set with Gaussian distribution, use "Gaussian"; For a data set with Ising distribution, use "Ising". Default is "Gaussian".}
  \item{rule}{rules to combine matrices that encode the conditional dependence relationships between sets of two observed variables. Options are "AND" and "OR". Default is "AND".}
}

\details{
The corlatent method has two assumptions. Assumption 1 states that the \eqn{R} replicates are assumed to follow a one-lag vector autoregressive model, conditioned on the latent variables. 
Assumption 2 states that the latent variables are piecewise constant across replicates.
Based on these two assumptions, the method solve the following problem for \eqn{1 \le j \le p}.
\deqn{
\min_{\theta_{j,-j}, \alpha_j, \Delta_j} \{ -\frac{1}{nR}l(\theta_{j,-j}, \alpha_j, \Delta_j) + \lambda\|\theta_{j,-j}\|_1 + \beta\|\alpha_j\|_1 + \gamma\|(I_n \otimes C)\Delta_j\|_1 \},
}
where \eqn{l(\theta_{j,-j}, \alpha_j, \Delta_j)} is the log likelihood function, \eqn{\theta_{j,-j}} encodes the conditional dependence relationships between \eqn{j}th observed variable and the other observed variables, \eqn{\alpha_j} models the correlation among replicates, \eqn{\Delta_j} encodes the latent effect, \eqn{\lambda}, \eqn{\beta}, \eqn{\gamma} are the tuning parameters, \eqn{I_n} is an n-dimensional identity matrix and \eqn{C} is the discrete first derivative matrix where the \eqn{i}th and \eqn{(i+1)}th column of every ith row are -1 and 1, respectively.
This method aims at modeling exponential family graphical models with correlated replicates and latent variables.
}

\value{
  \item{omega}{a matrix that encodes the conditional dependence relationships between sets of two observed variables}
  \item{theta}{the adjacency matrix with 0 and 1 encoding conditional independence and dependence between sets of two observed variables, respectively}
  \item{penalties}{the penalty values}
}

\references{
Jin, Y., Ning, Y., and Tan, K. M. (2020), `Exponential Family Graphical Models with Correlated Replicates and Unmeasured Confounders', preprint available.
}

\examples{
# Gaussian distribution with "AND" rule
n <- 20
R <- 10
p <- 5
l <- 2
s <- 2
seed <- 1

data <- generate_Gaussian(n, R, p, l, s, sparsityA = 0.95, sparsityobserved = 0.9,
sparsitylatent = 0.2, lwb = 0.3, upb = 0.3, seed)$X

result <- corlatent(data, accuracy = 1e-6, n, R, p,lambda1 = 0.1, lambda2 = 0.1,
lambda3 = 1e+5,distribution = "Gaussian", rule = "AND")
}
