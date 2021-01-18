#' @docType package
#' @name lavaSearch2
#' @aliases lavaSearch2, lavaSearch2-package
#'
#' @title Tools for Model Specification in the Latent Variable Framework
#' @description
#' The package contains three main functionalities:
#' \itemize{
#' \item \code{\link{compare2}}, \code{\link{summary2}}: Wald tests/robust Wald tests/F-tests/robust F-tests with improved control of the type 1 error in small samples.
#' \item \code{\link{glht2}}: adjustment for multiple comparisons when doing inference for multiple latent variable models.
#' \item \code{\link{modelsearch2}}: searching for local dependencies with adjustment for multiple comparisons.
#' }
#' It contains other useful functions such as:
#'\itemize{
#' \item \code{\link{calibrateType1}}: simulation study of the type 1 error of Wald tests.
#' \item \code{\link{createContrast}}: user-friendly function generating a contrast matrix.
#' \item \code{\link{getVarCov2}}: reconstruct the conditional variance covariance matrix.
#' \item \code{\link{iidJack}}: extract the jackknife iid decomposition.
#' }
#' @details
#' The latent variable models (LVM) considered in this package can be written \cr
#' as a measurement model:
#' \deqn{Y_i = \nu + \eta_i \Lambda + X_i K + \epsilon_i}
#' and a structural model:
#' \deqn{\eta_i = \alpha + \eta_i B + X_i \Gamma + \zeta_i}
#' where \eqn{\Sigma} is the variance covariance matrix of the residuals \eqn{\epsilon}, \cr
#' and \eqn{\Psi}   is the variance covariance matrix of the residuals \eqn{\zeta}. \cr
#' 
#' The corresponding conditional mean is:
#' \deqn{
#' \mu_i(\theta) = E[Y_i|X_i] = \nu + (\alpha + X_i \Gamma) (1-B)^{-1} \Lambda + X_i K
#' }
#' \deqn{
#' \Omega(\theta) = Var[Y_i|X_i] = \Lambda^{t} (1-B)^{-t} \Psi (1-B)^{-1} \Lambda + \Sigma
#' }
#'
#' The package aims to provides tool for testing linear hypotheses on the model coefficients
#' \eqn{\nu}, \eqn{\Lambda}, \eqn{K}, \eqn{\Sigma},
#' \eqn{\alpha}, \eqn{B}, \eqn{\Gamma}, \eqn{\Psi}.
#' Searching for local dependency enable to test whether the proposed model is too simplistic and if so to identify which additional coefficients should be added to the model.
#' 
#'
#'
#' 
#' @section Limitations:
#' 
#' 'lavaSearch2' has been design for Gaussian latent variable models.
#' This means that it may not work / give valid results:
#' \itemize{
#' \item in presence of censored or binary outcomes.
#' \item with stratified models (i.e. object of class \code{multigroup}).
#' } 
#' 

#' @useDynLib lavaSearch2, .registration=TRUE
#' @import lava
#' @importFrom ggplot2  aes_string autoplot
#' @importFrom graphics par plot text
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag
#' @importFrom methods is
#' @importFrom multcomp glht 
#' @importFrom mvtnorm pmvnorm qmvnorm rmvnorm qmvt pmvt
#' @importFrom parallel detectCores makeCluster stopCluster
#' @import Rcpp
#' @importFrom reshape2 melt
#' @importFrom sandwich estfun
#' @importFrom stats anova as.formula coef cov df.residual dist formula hclust logLik median model.frame model.matrix na.omit optim p.adjust pf pnorm predict qqnorm quantile pt residuals rnorm sd setNames sigma update vcov
#' @importFrom utils methods packageVersion setTxtProgressBar tail txtProgressBar
#' 
NULL



  
