#' @title bdpnormal Object Print
#' @description \code{print} method for class \code{bdpnormal}.
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x object of class \code{bdpnormal}. The result of a call to the
#'   \code{\link{bdpnormal}} function.
#'
#' @details Returns same output as a call to \code{\link[=summary,bdpnormal-method]{summary}}.
#' @export
setMethod("print", signature(x = "bdpnormal"), function(x){
  ### Return summary
  summary(x)
})


#' @title bdpbinomial Object Print
#' @description \code{print} method for class \code{bdpbinomial}.
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats sd density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x object of class \code{bdpbinomial}. The result of a call to the
#'   \code{\link{bdpbinomial}} function.
#'
#' @details Returns same output as a call to
#' \code{\link[=summary,bdpbinomial-method]{summary}}.
#' @export
setMethod("print", signature(x = "bdpbinomial"), function(x){
  ### Return summary
  summary(x)
})


#' @title bdpsurvival Object Print
#' @description \code{print} method for class \code{bdpsurvival}.
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x object of class \code{bdpsurvival}. The result of a call to the
#'   \code{\link{bdpsurvival}} function.
#' @details Displays a print of the \code{bdpsurvival} fit. The output
#'   is different, conditional on a one- or two-arm survival analysis.
#'
#'   In the case of a one-arm analysis, a brief summary is displayed,
#'   including the current data sample size, number of events,
#'   user input survival time, the augmented median survival probability,
#'   and corresponding lower and upper 95 percent interval limits.
#'
#'   When a control arm is present, the output is the same as a call to
#'   \code{\link[=summary,bdpsurvival-method]{summary}}.
#'
#' @export
setMethod("print", signature(x = "bdpsurvival"), function(x){
  posterior_treatment <- x$posterior_treatment
  posterior_control   <- x$posterior_control
  surv_time           <- x$args1$surv_time

  args1               <- x$args1
  data                <- args1$data
  data_current        <- args1$data_current
  breaks              <- args1$breaks
  arm2                <- args1$arm2


  treatment = NULL
  historical = NULL

  if(!arm2){
    ##############################################################################
    # Survival probability and surv_time
    ##############################################################################
    ### Print the augmented posterior
    survival_time_posterior_flat <- ppexp(surv_time,
                                          posterior_treatment$posterior_hazard,
                                          cuts=c(0,breaks))

    data_t <- subset(data, historical==0 & treatment == 1)
    n      <- nrow(data_t)
    s_t    <- with(data_t, Surv(time, status))# , type="mstate"))
    s_t    <- survival::survfitKM(factor(rep(1,n)), s_t)

    print_1arm <- matrix(c(nrow(data_current),
                         sum(s_t$n.event),
                         surv_time,
                         1-median(survival_time_posterior_flat),
                         1-quantile(survival_time_posterior_flat,0.975),
                         1-quantile(survival_time_posterior_flat,0.025)),nrow=1)
    print_1arm <- round(print_1arm,4)
    cnames <- c("n","events","surv_time","median","lower 95% CI","upper 95% CI")
    dimnames(print_1arm) <- list(rep("", nrow(print_1arm)), cnames)
    cat("\n")
    cat("    One-armed bdp survival\n\n")
    cat("\n")
    print(print_1arm)
  } else{
    ### Return summary
    summary(x)
  }

})




#' @title bdplm Object Print
#' @description \code{print} method for class \code{bdplm}.
#' @import methods
#' @importFrom utils head
#' @importFrom utils write.table
#' @importFrom stats density is.empty.model median model.offset model.response pweibull quantile rbeta rgamma rnorm var vcov
#' @param x object of class \code{bdplm}. The result of a call to the
#'   \code{\link{bdplm}} function.
#' @details Displays a print of the \code{bdplm} fit and the initial function call.
#'   The fit includes the estimate of the intercept, treatment effect, and
#'   covariate effects. The discount function weight estimates are displayed as well.
#'   If \code{method}="mc", the median estimate of alpha is displayed.
#'
#' @export
setMethod("print", signature(x = "bdplm"), function(x){

  # Format coefficients
  coefs <- x$estimates$coefficients
  p     <- ncol(coefs)
  coefs <- coefs[,-p]
  names(coefs)[1] <- "(Intercept)"
  coefs[1,] <- round(coefs[1,], 3)
  dimnames(coefs) <- list("", names(coefs))

  # Format alpha
  alpha_mat <- apply(x$alpha_discount, 2, median)
  alpha_mat <- matrix(alpha_mat, nrow=1)
  dimnames(alpha_mat) <- list("", names(x$alpha_discount))

  # Print output
  cat("\n")
  cat("Call:\n")
  print(x$args1$call)
  cat("\n\n")
  cat("Coefficients:\n")
  print(coefs)
  cat("\n\n")
  cat("Discount function value (alpha):\n")
  print(alpha_mat)
  cat("\n")
})
