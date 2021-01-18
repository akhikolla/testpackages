########################    Head     ########################

#' Print the Top Elements of the Joint CDF, Marginal Survival, and Conditional CDF after Nonparametric Analysis
#'
#' @description This function prints the top elements of the joint cdf, marginal survival, and conditional cdf from a \verb{bivrecNP} object.
#' @param x A \verb{bivrecNP} object
#' @param ... additional parameters if needed
#'
#' @export

head.bivrecNP <- function(x, ...) {

  if (!inherits(x, "bivrecNP")) stop("Must be a bivrecNP")

  cat("\nJoint CDF:\n", " ", sep = "")

  print(x$joint_cdf[1:6, ])

  cat("\nMarginal Survival:\n", " ", sep = "")

  print(x$marginal_survival[1:6, ])

  if (x$conditional==TRUE) {

    cat("\nConditional CDF:\n", " ", sep = "")

    print(x$conditional_cdf$conditional[1:6, ])
  }

}
