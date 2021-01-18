#' @template test_stat
#' @templateVar name Kolmogorov--Smirnov test statistic
#' @templateVar cite Miao2006
#' @templateVar formula \deqn{ \sup_t\left|F_n(t)-(1-F_n(-t))\right| }
#' @export
KS <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  KS_Cpp(X)
}

#' @template test_stat
#' @templateVar name The Sign test statistic
#' @templateVar cite Miao2006
#' @templateVar formula \deqn{ \frac1n \sum_{i=1}^nI\{X_i > 0\} - \frac12}
#' @export
SGN <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  SGN_Cpp(X)
}

#' @template test_stat
#' @templateVar name The Wilcoxon test statistic
#' @templateVar cite Miao2006
#' @templateVar formula \deqn{ \frac{1}{{n\choose2}} \sum_{1\leq i<j\leq n} I\{X_i+X_j> 0\} - \frac12 }
#' @export
WCX <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  WCX_Cpp(X)
}
