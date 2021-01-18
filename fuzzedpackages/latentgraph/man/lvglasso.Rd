\name{lvglasso}
\alias{lvglasso}

\title{
Estimate Gaussian Graphical Models with Latent Variables
}
\description{
Estimate Gaussian graphical models with latent variables using the method in Chandrasekaran et al. (2012).
}
\usage{
lvglasso(data, n, p, lambda1, lambda2, rule = "AND")
}

\arguments{
  \item{data}{data set, can be a matrix or data frame with \eqn{n} rows and \eqn{p} columns.}
  \item{n}{the number of observations.}
  \item{p}{the number of observed variables.}
  \item{lambda1}{tuning parameter that encourages estimated graph to be sparse. Lambda1 is proportional to lambda2.}
  \item{lambda2}{tuning parameter that encourages the matrix \eqn{K_{O,H} (K_H)^{-1} K_{H, O}} to be low rank, where \eqn{K_H} and \eqn{K_{O,H}} quantify the dependecies among the latent variables and between the observed variables and latent variables, respectively. The matrix \eqn{K_{O,H} (K_H)^{-1} K_{H, O}} summarizes the impact of marginalization over the latent variables.}
  \item{rule}{rules to combine the inverse covariance matrices. Options are "AND" and "OR". Default is "AND".}
}

\details{
The lvglasso method assumes that all the variables, both observed and latent, are jointly Gaussian, and specifies the conditional distribution of observed variables on the latent variables by a graphical model. Under the high-dimentional setting, this method provides consistent estimators for the conditional graphical model of observed variables conditioned on latent variables. 
}

\value{
  \item{omega}{a matrix that encodes the conditional dependence relationships between sets of two observed variables}
  \item{theta}{the adjacency matrix with 0 and 1 encoding conditional independence and dependence between sets of two observed variables, respectively}
  \item{penalties}{the penalty values}
}

\references{
Chandrasekaran, V., Parrilo, P. A. & Willsky, A. S. (2012), `Latent variable graphical model selection via convex optimization', Ann. Statist. 40(4), 1935--1967.
}


\examples{
#Gaussian distribution with "AND" rule
n <- 50
R <- 20
p <- 30
l <- 2
s <- 2
data <- generate_Gaussian(n, R, p, l, s, sparsityA = 0.95, sparsityobserved = 0.9,
sparsitylatent = 0.2, lwb = 0.3, upb = 0.3, seed = 1)$X

result <- lvglasso(data, n, p, lambda1 = 0.222, lambda2 = 0.1*0.222, rule = "AND")
}

