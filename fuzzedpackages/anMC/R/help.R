#' @title anMC package
#' @description Efficient estimation of high dimensional orthant probabilities. The package main functions are: \itemize{
#'     \item \code{\link{ProbaMax}}: the main function for high dimensional othant probabilities. Computes \eqn{P(max X > t)}, where \eqn{X} is a Gaussian vector and \eqn{t} is the selected threshold. It implements the \code{GANMC} algorithm and allows for user-defined sampler and core probability estimates.
#'     \item \code{\link{ProbaMin}}: analogous of \code{ProbaMax} for the problem \eqn{P(min X < t)}, where \eqn{X} is a Gaussian vector and \eqn{t} is the selected threshold.  It implements the \code{GANMC} algorithm and allows for user-defined sampler and core probability estimates.
#'     \item \code{\link{conservativeEstimate}}: the main function for conservative estimates computation. Requires the mean and covariance of the posterior field at a discretization design.
#' }
#' @details Package: anMC \cr
#' Type: Package \cr
#' Version: 0.2.2 \cr
#' Date: 2019-10-23
#'
#' @author Dario Azzimonti (dario.azzimonti@@gmail.com) . Thanks to David Ginsbourger for the fruitful discussions and his continuous help in testing and improving the package.
#' @docType package
#' @name anMC
#' @import mvtnorm
#' @importFrom Rcpp evalCpp
#' @importFrom stats cov dist lm pnorm quantile var
#' @useDynLib anMC
#' @note This work was supported in part by the Swiss National Science Foundation, grant number 146354 and the Hasler Foundation, grant number 16065.
#' @references Azzimonti, D. and Ginsbourger, D. (2018). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Journal of Computational and Graphical Statistics, 27(2), 255-267. \href{https://doi.org/10.1080/10618600.2017.1360781}{DOI: 10.1080/10618600.2017.1360781}
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#'
#' Bolin, D. and Lindgren, F. (2015). Excursion and contour uncertainty regions for latent Gaussian models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(1):85--106.
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Dickmann, F. and Schweizer, N. (2014). Faster comparison of stopping times by nested conditional Monte Carlo. arXiv preprint arXiv:1402.0243.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#'
#' Genz, A. and Bretz, F. (2009). Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics 195. Springer-Verlag.
#'
#' Horrace, W. C. (2005). Some results on the multivariate truncated normal distribution. Journal of Multivariate Analysis, 94(1):209--221.
#'
#' Robert, C. P. (1995). Simulation of truncated normal variables. Statistics and Computing, 5(2):121--125.
NULL
