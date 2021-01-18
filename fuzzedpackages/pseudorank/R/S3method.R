################################################################################
### File: S3method.R
### Description: S3 generic functions to handle vectors and data.frames / formula objects
###              for psrank and hettmansperger_norton_test
###
################################################################################



#' Calculation of Pseudo-Ranks
#'
#' @description Calculation of (mid) pseudo-ranks of a sample. In case of ties (i.e. equal values), the average of min pseudo-ranks and max-pseudo-ranks are taken (similar to rank with ties.method="average").
#' @param x vector containing the observations
#' @param y vector specifiying the group to which the observations from the x vector belong to
#' @param data data.frame containing the variables in the formula (observations and group)
#' @param formula formula object
#' @param na.last for controlling the treatment of NAs. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed (recommended).
#' @param ties.method type of pseudo-ranks: either 'average' (recommended), 'min' or 'max'.
#' @param ... further arguments
#' @return Returns a numerical vector containing the pseudo-ranks.
#' @rdname pseudorank
#' @references Brunner, E., Bathke, A.C., and Konietschke, F. (2018a). Rank- and Pseudo-Rank Procedures for Independent Observations in Factorial Designs - Using R and SAS. Springer Series in Statistics, Springer, Heidelberg. ISBN: 978-3-030-02912-8.
#' @references Happ M, Zimmermann G, Brunner E, Bathke AC (2020). Pseudo-Ranks: How to Calculate Them Efficiently in R. Journal of Statistical Software, Code Snippets, *95*(1), 1-22. doi: 10.18637/jss.v095.c01 (URL:https://doi.org/10.18637/jss.v095.c01).
#' @example R/example_1.txt
#' @keywords export
pseudorank <- function(x, ...){
  UseMethod("pseudorank")
}

#' @method pseudorank numeric
#' @rdname pseudorank
#' @keywords export
pseudorank.numeric <- function(x, y, na.last = NA, ties.method = c("average", "max", "min"), ...){
  stopifnot(na.last %in% c(TRUE, FALSE, NA))
  ties.method = match.arg(ties.method)
  recursiveCalculation(x, y, na.last, ties.method)
}

#' @method pseudorank formula
#' @rdname pseudorank
#' @keywords export
pseudorank.formula <- function(formula, data, na.last = NA, ties.method = c("average", "max", "min"), ...){
  stopifnot(na.last %in% c(TRUE, FALSE, NA))
  ties.method = match.arg(ties.method)
  df <- model.frame(formula, data, na.action = NULL)
  recursiveCalculation(df[, 1], df[, 2], na.last, ties.method)
}

#' Calculation of Pseudo-Ranks (Deprecated)
#'
#' @description Calculation of (mid) pseudo-ranks of a sample. In case of ties (i.e. equal values), the average of min pseudo-ranks and max-pseudo-ranks are taken (similar to rank with ties.method="average").
#' @param x vector containing the observations
#' @param ... further arguments (see help for pseudorank)
#' @return Returns a numerical vector containing the pseudo-ranks.
#' @rdname psrank-deprecated
#' @example R/example_1.txt
#' @references Happ M, Zimmermann G, Brunner E, Bathke AC (2020). Pseudo-Ranks: How to Calculate Them Efficiently in R. Journal of Statistical Software, Code Snippets, *95*(1), 1-22. doi: 10.18637/jss.v095.c01 (URL:https://doi.org/10.18637/jss.v095.c01).
#' @keywords export
psrank <- function(x, ...) {
  .Deprecated("pseudorank", package=NULL,
              old = as.character(sys.call(sys.parent()))[1L])
  UseMethod("pseudorank")
}


#' Hettmansperger-Norton Trend Test for k-Samples
#'
#' @description This function calculates the Hettmansperger-Norton trend test using pseudo-ranks under the null hypothesis H0F: F_1 = ... F_k.
#' @rdname hettmansperger_norton_test
#' @param x vector containing the observations
#' @param y vector specifiying the group to which the observations from the x vector belong to
#' @param data data.frame containing the variables in the formula (observations and group)
#' @param formula formula object
#' @param na.rm a logical value indicating if NA values should be removed
#' @param alternative either decreasing (trend k, k-1, ..., 1) or increasing (1, 2, ..., k) or custom (then argument trend must be specified)
#' @param trend custom numeric vector indicating the trend for the custom alternative, only used if alternative = "custom"
#' @param pseudoranks logical value indicating if pseudo-ranks or ranks should be used
#' @param ... further arguments are ignored
#' @return Returns an object.
#' @example R/example_2.txt
#' @references Brunner, E., Bathke, A.C., and Konietschke, F. (2018a). Rank- and Pseudo-Rank Procedures for Independent Observations in Factorial Designs - Using R and SAS. Springer Series in Statistics, Springer, Heidelberg. ISBN: 978-3-030-02912-8.
#' @references Happ M, Zimmermann G, Brunner E, Bathke AC (2020). Pseudo-Ranks: How to Calculate Them Efficiently in R. Journal of Statistical Software, Code Snippets, *95*(1), 1-22. doi: 10.18637/jss.v095.c01 (URL:https://doi.org/10.18637/jss.v095.c01).
#' @references Hettmansperger, T. P., & Norton, R. M. (1987). Tests for patterned alternatives in k-sample problems. Journal of the American Statistical Association, 82(397), 292-299
#' @keywords export
hettmansperger_norton_test <- function(x, ...) {
  UseMethod("hettmansperger_norton_test")
}

#' @method hettmansperger_norton_test numeric
#' @rdname hettmansperger_norton_test
#' @keywords export
hettmansperger_norton_test.numeric <- function(x, y, na.rm = FALSE, alternative = c("decreasing", "increasing", "custom"), trend = NULL, pseudoranks = TRUE, ...) {
  alternative = match.arg(alternative)
  return(hettmansperger_norton_test_internal(x, y, na.rm, alternative = alternative, formula = NULL, trend = trend, pseudoranks = pseudoranks, ...))
}

#' @method hettmansperger_norton_test formula
#' @rdname hettmansperger_norton_test
#' @keywords export
hettmansperger_norton_test.formula <- function(formula, data, na.rm = FALSE, alternative = c("decreasing", "increasing", "custom"), trend = NULL, pseudoranks = TRUE, ...) {
  alternative = match.arg(alternative)
  model <- model.frame(formula, data = data, na.action = NULL)
  colnames(model) <- c("data", "group")
  return(hettmansperger_norton_test_internal(model$data, model$group, na.rm, alternative = alternative, formula = formula, trend = trend, pseudoranks = pseudoranks, ...))
}







#' Kruskal-Wallis Test
#'
#' @description This function calculates the Kruskal-Wallis test using pseudo-ranks under the null hypothesis H0F: F_1 = ... F_k.
#' @rdname kruskal_wallis_test
#' @param x numeric vector containing the data
#' @param grp factor specifying the groups
#' @param na.rm a logical value indicating if NA values should be removed
#' @param formula optional formula object
#' @param data optional data.frame of the data
#' @param pseudoranks logical value indicating if pseudo-ranks or ranks should be used
#' @param ... further arguments are ignored
#' @return Returns an object of class 'pseudorank'
#' @example R/example_3.txt
#' @references Brunner, E., Bathke, A.C., and Konietschke, F. (2018a). Rank- and Pseudo-Rank Procedures for Independent Observations in Factorial Designs - Using R and SAS. Springer Series in Statistics, Springer, Heidelberg. ISBN: 978-3-030-02912-8.
#' @keywords export
kruskal_wallis_test <- function(x, ...) {
  UseMethod("kruskal_wallis_test")
}

#' @method kruskal_wallis_test numeric
#' @rdname kruskal_wallis_test
#' @keywords export
kruskal_wallis_test.numeric <- function(x, grp, na.rm = FALSE, pseudoranks = TRUE, ...) {
  return(kruskal_wallis_internal(data=x, group=grp, na.rm = na.rm, pseudoranks = pseudoranks, formula = NULL, ...))
}

#' @method kruskal_wallis_test formula
#' @rdname kruskal_wallis_test
#' @keywords export
kruskal_wallis_test.formula <- function(formula, data, na.rm = FALSE, pseudoranks = TRUE, ...) {
  model <- model.frame(formula, data = data, na.action = NULL)
  colnames(model) <- c("data", "group")
  return(kruskal_wallis_internal(model$data, as.factor(model$group), na.rm = na.rm, formula = formula, pseudoranks = pseudoranks, ...))
}





#' Kepner-Robinson Test
#'
#' @description This function calculates the Kepner-Robinson test using ranks under the null hypothesis H0F: F_1 = ... F_k where F_i are the marginal distributions. Each subject needs to have k measurements. This test assumes that the covariance matrix of a subject has a compound symmetry structure.
#' @rdname kepner_robinson_test
#' @param x numeric vector containing the data
#' @param time factor specifying the groups
#' @param subject factor specifying the subjects or the name of the subject column if a data.frame is used
#' @param na.rm a logical value indicating if NA values should be removed
#' @param formula optional formula object
#' @param data optional data.frame of the data
#' @param distribution either 'Chisq' or 'F' approximation
#' @param ... further arguments are ignored
#' @return Returns an object of class 'pseudorank'
#' @example R/example_3.txt
#' @references James L. Kepner & David H. Robinson (1988) Nonparametric Methods for Detecting Treatment Effects in Repeated-Measures Designs, Journal of the American Statistical Association, 83:402, 456-461.
#' @keywords internal
kepner_robinson_test <- function(x, ...) {
  UseMethod("kepner_robinson_test")
}

#' @method kepner_robinson_test numeric
#' @rdname kepner_robinson_test
#' @keywords internal
kepner_robinson_test.numeric <- function(x, time, subject, na.rm = FALSE, distribution = c("Chisq", "F"), ...) {
  distribution = match.arg(distribution)
  return(kepner_robinson_test_internal(data=x, time=as.factor(time), subject=as.factor(subject), na.rm = na.rm, distribution = distribution, formula = NULL, ...))
}

#' @method kepner_robinson_test formula
#' @rdname kepner_robinson_test
#' @keywords internal
kepner_robinson_test.formula <- function(formula, data, subject, na.rm = FALSE, distribution = c("Chisq", "F"), ...) {
  stopifnot(is.character(subject))
  distribution = match.arg(distribution)
  model <- model.frame(formula, data = data, na.action = NULL)
  model$subject <- data[, subject]
  colnames(model) <- c("data", "time", "subject")
  return(kepner_robinson_test_internal(data=model$data, time=as.factor(model$time), subject=as.factor(model$subject), na.rm = na.rm, distribution = distribution, formula = formula, ...))
}





#' @keywords export
print.pseudorank <- function(x, ...) {
  cat(x$name)
  cat("\n","\n")
  if(!is.null(x$formula)) {
    cat("Call:", "\n")
    print(x$formula)
    cat("\n")
  }
  if(!is.null(x$alternative)){
    cat("Alternative: ", x$alternative, "\n")
    if(x$alternative == "custom") {
      cat("Trend: ", x$trend, "\n")
    }
  }
  cat("Test Statistic: ", x$test, "\n")
  cat("Distribution of Statistic: ", x$distribution, "\n")
  if(!is.null(x$df)) {
    cat("Degrees of Freedom: ", paste(x$df,collapse=", "), "\n")
  }
  cat("unweighted relative Effects / Pseudo-ranks: ", x$pseudoranks)
  cat("\n")
  cat("p-Value: ", x$pValue, "\n")
  # cat("\n")
  # cat("Descriptive:\n")
  # df <- data.frame(n = x$ss, p = x$pHat)
  # print(df, row.names = FALSE)
}

#' @keywords export
summary.pseudorank <- function(object, ...) {
  cat(object$name)
  cat("\n", "\n")
  if(!is.null(object$formula)) {
    cat("Call:", "\n")
    print(object$formula)
    cat("\n")
  }
  if(!is.null(object$alternative)){
    cat("Alternative: ", object$alternative, "\n")
    if(object$alternative == "custom") {
      cat("Trend: ", object$trend, "\n")
    }
  }
  cat("Test Statistic: ", object$test, "\n")
  cat("Distribution of Statistic: ", object$distribution, "\n")
  if(!is.null(object$df)) {
    cat("Degrees of Freedom: ", paste(object$df,collapse=", "), "\n")
  }
  cat("unweighted relative Effects / Pseudo-ranks: ", object$pseudoranks)
  cat("\n")
  cat("p-Value: ", object$pValue, "\n")
  cat("\n")
  cat("Descriptive:\n")
  df <- data.frame(n = object$ss, p = object$pHat)
  print(df, row.names = FALSE)
}

