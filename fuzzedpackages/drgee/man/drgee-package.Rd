\encoding{latin1}
\name{drgee-package}
\alias{drgee-package}
\docType{package}
\title{
Doubly Robust Generalized Estimating Equations
}
\description{
  The main function is \code{drgee}, which estimates a parameter
  \eqn{\beta}{beta} in a model defined as
  \eqn{g\{E(Y|A,L)\}-g\{E(Y|A=0,L)\}=\beta^T \{A\cdot X(L)\}}{g\{E(Y|A,L)\}-g\{E(Y|A=0,L)\}=beta^T
    (A * X(L))\}}. By supplying (nuisance) models for
  \eqn{g\{E(Y|A=0,L)}
  and \eqn{E(A|L)}, a consistent estimate of \eqn{\beta}{beta} is
  obtained when at least one of these models is correctly specified.

  In the function \code{drgee}, three estimation methods are
  implemented:
  O-estimation, where a model for \eqn{g\{E(Y|A=0,L)} is used;
  E-estimation, where a model for \eqn{E(A|L)} is used; and
  DR-estimation where models for both \eqn{g\{E(Y|A=0,L)} and
  \eqn{E(A|L)} are used.

  The function \code{gee} is an implementation of standard GEE with
  independent working correlation matrix.

  For conditional methods with clustered data, cluster-specific
  intercepts are assumed in all models.

  The function \code{drgeeData} is used for extraction and manipulation
  of data, and is called by \code{drgee} and \code{gee}.

  The function \code{RobVcov} is used to calculate standard errors of
  the estimates, given a vector of the residuals from estimating
  equations, the Jacobian, and a cluster-identifying variable.

  The function \code{findRoots} solves a system of non linear equations.
}
\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander} with contributions from Alexander Ploner.
}
\keyword{models}
\keyword{regression}
\keyword{cluster}
