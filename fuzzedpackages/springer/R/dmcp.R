#' The first order derivative function of MCP (Minimax Concave Penalty)
#' @param theta a coefficient vector.
#' @param lambda the tuning parameter.
#' @param gamma the regularization parameter for MCP (Minimax Concave Penalty).
#' It balances the unbiasedness and concavity of MCP.
#' @return the first order derivative of the MCP function.
#' @details
#' The regularization parameter \eqn{\gamma} for MCP should be obtained via a data-driven approach in a rigorous way.
#' Among the published studies, it is suggested to check several choices, such as 1.4, 3, 4.2, 5.8, 6.9, and 10, then fix the value.
#' We examined this sequence in our study and found that the results are not sensitive to the choice of value for \eqn{\gamma}.  Therefore,
#' we set the value to 3. To be prudent, other values should also be examined in practice.
#'
#' @references
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang, Y. and Wu, C. (2019). Robust network-based regularization and variable selection for high-dimensional genomic data in cancer prognosis.
#' \href{https://doi.org/10.1002/gepi.22194}{\emph{Genetic epidemiology}, 43(3), 276-291}
#'
#' Wu, C., Zhang, Q., Jiang, Y. and Ma, S. (2018). Robust network-based analysis of the associations between (epi) genetic measurements.
#' \href{https://doi.org/10.1016/j.jmva.2018.06.009}{\emph{Journal of multivariate analysis}, 168, 119-130}
#'
#' Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y. and Wu, C. (2017). Network-based regularization for high dimensional SNP data in the case–control study of Type 2 diabetes.
#' \href{https://doi.org/10.1186/s12863-017-0495-5}{\emph{BMC genetics}, 18(1), 44}
#'
#' @examples
#' theta=runif(30,-4,4)
#' lambda=1
#' dmcp(theta,lambda,gamma=3)
#' @export

dmcp <- function(theta,lambda,gamma=3){
  #length of parameter
  p<-length(theta)
  #create vector of zeros
  b1<-rep(0,p)
  #if theta is less than gamma*lambda set it to 1
  b1[abs(theta)<=(gamma*lambda)]<-1
  (lambda-theta/gamma)*b1
}
