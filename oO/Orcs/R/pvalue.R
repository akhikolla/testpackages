#' Get \emph{p}-value from 'lm' object
#' 
#' @description 
#' Retrieve the \emph{p}-value associated with a univariate linear regression. 
#' 
#' @param mod An object of class \code{lm}.
#' 
#' @return 
#' A \code{numeric} \emph{p}-value.
#' 
#' @source  
#' \href{https://stat.ethz.ch/pipermail/r-help/2005-December/084308.html}{retrieving p-values in lm}
#' on R-help mailing list (last accessed 2018-06-02).
#' 
#' @seealso 
#' \code{\link{lm}}.
#' 
#' @examples 
#' ## taken from ?lm
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' 
#' pvalue(lm.D9)
#' 
#' @export pvalue
#' @name pvalue
pvalue <- function(mod) {
  
  mod_sum <- summary(mod)
  coef(mod_sum)[2, 4]
  
}