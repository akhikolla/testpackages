##### CLASS DEFINITIONS ########
#' @docType class
#' @title cokm Class
#' @description This is an S4 class definition for \code{\link{cokm}} in the \code{\link{ARCokrig}} package
#' 
#' @slot output a list of \eqn{s} elements, each of which contains a matrix of computer model outputs.
#' @slot input a list of \eqn{s} elements, each of which contains a matrix of inputs.
#' @slot param a list of \eqn{s} elements, each of which contains a vector of initial values for 
#' correlation parameters (and nugget variance parameters if 
#' nugget terms are included in AR-cokriging models).
#' @slot cov.model a string indicating the type of covariance
#' function in AR-cokriging models. Current covariance functions include
#' \describe{
#' \item{exp}{product form of exponential covariance functions.} 
#' \item{matern_3_2}{product form of Matern covariance functions with 
#' smoothness parameter 3/2.}
#' \item{matern_5_2}{product form of Matern covariance functions with
#' smoothness parameter 5/2.}
#' \item{Gaussian}{product form of Gaussian covariance functions.}
#' \item{powexp}{product form of power-exponential covariance functions with roughness parameter fixed at 1.9.}
#' }
#' 
#' @slot nugget.est a logical value indicating whether nugget parameter is included or not. Default value is \code{FALSE}.
#' @slot prior a list of arguments to setup the prior distributions with the reference prior as default
#' \describe{
#'	\item{name}{the name of the prior. Current implementation includes 
#'  \code{JR}, \code{Reference}, \code{Jeffreys}, \code{Ind_Jeffreys}}
#'  \item{hyperparam}{hyperparameters in the priors. 
#'	For jointly robust (JR) prior, three parameters are included: 
#' \eqn{a} refers to the polynomial penalty to avoid singular correlation 
#'   matrix with a default value 0.2; \eqn{b} refers to the exponenetial penalty to avoid 
#'   diagonal correlation matrix with a default value 1; nugget.UB is the upper
#' bound of the nugget variance with default value 1, which indicates that the
#' nugget variance has support (0, 1).}
#'
#'}
#' 
#' @slot opt a list of arguments to setup the \code{\link{optim}} routine.
#' @slot NestDesign a logical value indicating whether the 
#' experimental design is hierarchically nested within each level
#' of the code.
#' 
#' @slot tuning a list of arguments to control the MCEM algorithm for non-nested
#' design. It includes the arguments 
#' \describe{
#'    \item{maxit}{the maximum number of MCEM iterations.}
#'    \item{tol}{a tolerance to stop the MCEM algorithm. If the parameter 
#'    difference between any two consecutive MCEM algorithm is less than 
#'    this tolerance, the MCEM algorithm is stopped.}
#'    \item{n.sample}{the number of Monte Carlo samples in the 
#'    MCEM algorithm.}
#'    \item{verbose}{a logical value to show the MCEM iterations if it is true.}
#'}
#'
#'
#' @slot info a list that contains 
#' \describe{
#'  \item{iter}{number of iterations used in the MCEM algorithm}
#'  \item{eps}{parameter difference after the MCEM algorithm stops.}
#'}
#' @keywords AR-Cokriging Objective-Bayes Computer-Experiments Uncertainty-Quantification
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
setClass("cokm", representation(
  formula = "list",
  output = "list",
  input = "list",
  param = "list",
  cov.model = "character",
  nugget.est = "logical",
  prior = "list",
  opt = "list",
  NestDesign = "logical",
  tuning = "list",
  info = "list")
) 





##### CLASS DEFINITIONS ########
#' @docType class
#' @title mvcokm Class
#' @description This is an S4 class definition for \code{\link{mvcokm}} in the \code{\link{ARCokrig}} package
#' 
#' @slot output a list of \eqn{s} elements, each of which contains a matrix of computer model outputs.
#' @slot input a list of \eqn{s} elements, each of which contains a matrix of inputs.
#' @slot param a list of \eqn{s} elements, each of which contains a vector of initial values for 
#' correlation parameters (and nugget variance parameters if 
#' nugget terms are included in AR-cokriging models).
#' @slot cov.model a string indicating the type of covariance
#' function in AR-cokriging models. Current covariance functions include
#' \describe{
#' \item{exp}{product form of exponential covariance functions.} 
#' \item{matern_3_2}{product form of Matern covariance functions with 
#' smoothness parameter 3/2.}
#' \item{matern_5_2}{product form of Matern covariance functions with
#' smoothness parameter 5/2.}
#' \item{Gaussian}{product form of Gaussian covariance functions.}
#' \item{powexp}{product form of power-exponential covariance functions with roughness parameter fixed at 1.9.}
#' \item{aniso_exp}{anisotropic form of exponential covariance function.}
#' \item{aniso_matern_3_2}{anisotropic form of Matern covariance functions with 
#' smoothness parameter 3/2.}
#' \item{aniso_matern_5_2}{anisotropic form of Matern covariance functions with 
#' smoothness parameter 5/2.}
#' }
#' 
#' @slot nugget.est a logical value indicating whether the nugget is included or not. Default value is \code{FALSE}.
#' 
#' @slot prior a list of arguments to setup the prior distributions with the jointly robust prior as default
#' \describe{
#'	\item{name}{the name of the prior. Current implementation includes 
#'  \code{JR}, \code{Reference}, \code{Jeffreys}, \code{Ind_Jeffreys}}
#'  \item{hyperparam}{hyperparameters in the priors. 
#'	For jointly robust (JR) prior, three parameters are included: 
#' \eqn{a} refers to the polynomial penalty to avoid singular correlation 
#'   matrix with a default value 0.2; \eqn{b} refers to the exponenetial penalty to avoid 
#'   diagonal correlation matrix with a default value 1; nugget.UB is the upper
#' bound of the nugget variance with default value 1, which indicates that the
#' nugget variance has support (0, 1).}
#'}
#' 
#' 
#' @slot opt a list of arguments to setup the \code{\link{optim}} routine.
#' @slot NestDesign a logical value indicating whether the 
#' experimental design is hierarchically nested within each level
#' of the code.
#' 
#' @slot tuning a list of arguments to control the MCEM algorithm for non-nested
#' design. It includes the arguments 
#' \describe{
#'    \item{maxit}{the maximum number of MCEM iterations.}
#'    \item{tol}{a tolerance to stop the MCEM algorithm. If the parameter 
#'    difference between any two consecutive MCEM algorithm is less than 
#'    this tolerance, the MCEM algorithm is stopped.}
#'    \item{n.sample}{the number of Monte Carlo samples in the 
#'    MCEM algorithm.}
#'}
#'
#' @slot info a list that contains 
#' \describe{
#'  \item{iter}{number of iterations used in the MCEM algorithm}
#'  \item{eps}{parameter difference after the MCEM algorithm stops}
#'}
#'
#' @keywords AR-Cokriging Objective-Bayes Computer-Experiments Uncertainty-Quantification
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
setClass("mvcokm", representation(
  formula = "list",
  output = "list",
  input = "list",
  param = "list",
  cov.model = "character",
  nugget.est = "logical",
  prior = "list",
  opt = "list",
  NestDesign = "logical",
  tuning="list",
  info="list")
) 


