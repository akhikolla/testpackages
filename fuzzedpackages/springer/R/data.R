#' simulated data for demonstrating the usage of springer
#'
#' Simulated gene expression data for demonstrating the usage of springer.
#'
#' @docType data
#' @keywords datasets
#' @name dat
#' @aliases dat
#' @usage data("dat")
#' @format The dat file consists of five components: e, g, y, clin and coeff. The coefficients are the true values of parameters used for generating Y.
#' @details
#'
#' \strong{The data model for generating Y}
#'
#' Consider a longitudinal case study with \eqn{n} subjects and \eqn{k_i} measurements over time for the \eqn{i}th subject (\eqn{i=1,\ldots,n}).
#' Let \eqn{Y_{ij}} be the response of the \eqn{j}th observation for the \eqn{i}th subject (\eqn{i=1,\ldots,n}, \eqn{j=1,\ldots,k_i}),
#' \eqn{X_{ij}=(X_{ij1},...,X_{ijp})^\top} be a \eqn{p}-dimensional vector of covariates denoting \eqn{p} genetic factors, \eqn{E_{ij}=(E_{ij1},...,E_{ijq})^\top}
#' be a \eqn{q}-dimensional environmental factor and \eqn{Clin_{ij}=(Clin_{ij1},...,Clin_{ijt})^\top} be a \eqn{t}-dimensional clinical factor.  There is time dependence among measurements on the same subject, but we assume that the measurements
#' between different subjects are independent.  The model we used for hierarchical variable selection for gene--environment interactions is given as:
#'
#' \deqn{Y_{ij}= \alpha_0 + \sum_{m=1}^{t}\theta_m Clin_{ijm} + \sum_{u=1}^{q}\alpha_u E_{iju} + \sum_{v=1}^{p}(\gamma_v X_{ijv} + \sum_{u=1}^{q}h_{uv} E_{iju} X_{ijv})+\epsilon_{ij},}
#' where \eqn{\alpha_{0}} is the intercept and the marginal density of \eqn{Y_{ij}} belongs to a canonical exponential family defined in Liang and Zeger (1986).
#' Define \eqn{\eta_v=(\gamma_v, h_{1v}, ..., h_{qv})^\top}, which is a vector of length q+1 and \eqn{Z_{ijv}=(X_{ijv}, E_{ij1}X_{ijv}, ..., E_{ijq}X_{ijv})^\top},
#' which contains the main genetic effect of the \eqn{v}th SNP from the \eqn{j}th measurement on the \eqn{i}th subject and its interactions with all the \eqn{q}
#' environmental factors. The model can be written as:
#'  \deqn{Y_{ij}= \alpha_0 + \sum_{m=1}^{t}\theta_m Clin_{ijm} + \sum_{u=1}^{q}\alpha_u E_{iju} + \sum_{v=1}^{p}\eta_v^\top Z_{ijv}+\epsilon_{ij},}
#' where \eqn{Z_{ijv}} is the \eqn{v}th genetic factor and its interactions with the \eqn{q} environment factors for the \eqn{j}th measurement on the \eqn{i}th subject,
#' and \eqn{\eta_{v}} is the corresponding coefficient vector of length \eqn{1+q}. The random error \eqn{\epsilon_{i}=(\epsilon_{i1},...,\epsilon_{ik_i})^{T}}, which is assumed to follow a multivariate normal distribution with \eqn{\Sigma_i}
#' as the covariance matrix for the repeated measurements of the \eqn{ith} subject among the \eqn{k_i} time points.
#'
#' @seealso \code{\link{springer}}
NULL
