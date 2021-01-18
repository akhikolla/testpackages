##########################################################################
#'
#' @title Test for conditional independence in the multivariate normal
#'     distribution
#' @description Test for conditional independence in the multivariate
#'     normal distribution.
#' @name citest-mvn
#'
##########################################################################
#'
#' @param x A list with elements \code{cov} and \code{n.obs} (such as returned
#'     from calling \code{cov.wt()} on a dataframe. See examples below.)
#' @param set A specification of the test to be made. The tests are of the form
#'     u and v are independent condionally on S where u and v are variables and
#'     S is a set of variables. See 'details' for details about specification of
#'     \code{set}.
#' @param statistic The test statistic to be used, valid choices are
#'     \code{"DEV"} and \code{"F"}.
#' @param \dots Additional arguments
#'
#' @details
#' 
#' \code{set} can be 1) a vector or 2) a right-hand sided formula in which
#' variables are separated by '+'. In either case, it is tested if the first
#' two variables in the \code{set} are conditionally independent given the
#' remaining variables in \code{set}.  (Notice an abuse of the '+' operator in
#' the right-hand sided formula: The order of the variables does matter.)
#' 
#' If \code{set} is \code{NULL} then it is tested whether the first two
#' variables are conditionally independent given the remaining variables.
#' 
#' \code{x} must be a list with components \code{cov} and \code{n.obs} such as
#' returned by calling \code{cov.wt( , method='ML')} on a dataframe.
#' 
#' @return An object of class `citest` (which is a list).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest}}, \code{\link{ciTest_table}},
#'     \code{\link{ciTest_df}},
#'     \code{\link{ciTest_mvn}}, \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#' 
#' data(carcass)
#' ciTest(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
#' ciTest_mvn(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
#' 
#' @export 
ciTest_mvn <- function(x, set=NULL, statistic="DEV", ...){

    if(any(is.na(match(c("cov", "n.obs"), names(x))))){
    stop("Expecting a list with components 'cov' and 'n.obs'\n")
  }

  if (is.null(set)){
    set   <- colnames(x$cov)
    x$cov <- x$cov[set, set]
  } else {
    if (inherits(set, c("formula", "character"))){
      set <- unlist(rhsFormula2list(set))
      set <- colnames(x$cov)[pmatch(set, colnames(x$cov))]
      x$cov <- x$cov[set, set]
    }
  }
  .ciTest_mvn_internal(x, statistic=statistic, ...) ## Should <set> go here...
}

### This is the workhorse.
###
.ciTest_mvn_internal <- function(x, statistic="DEV", ...){

  statistic <- match.arg(toupper(statistic), c("DEV", "F"))
  method <- if (identical(statistic, "DEV")) "CHISQ" else "F"

  S     <- x$cov
  n.obs <- x$n.obs

  vn <- colnames(S)
  K  <- length(vn)

  R   <- vn[-(1:2)]
  v1R <- c(vn[1],R)
  v2R <- c(vn[2],R)

  v1R.idx <- match(v1R, vn)
  v2R.idx <- match(v2R, vn)
  R.idx   <- match(R,   vn)

  d <- n.obs * (log(det(S[v1R.idx, v1R.idx, drop=FALSE])) +
                log(det(S[v2R.idx, v2R.idx, drop=FALSE])) -
                log(det(S[R.idx, R.idx, drop=FALSE])) - log(det(S))
                )

  num.df <- 1
  switch(statistic,
         "DEV"={
           tobs     <- d
           denom.df <- NULL
           p        <- 1 - pchisq(tobs, df=num.df)
         },
         "F"={
           tobs     <- (exp(d / n.obs) - 1) * (n.obs - K)
           denom.df <- n.obs - K
           p        <- 1 - pf(tobs, num.df, denom.df)
         })

  ans <- list(statistic=tobs, p.value=p, df=num.df, denom.df=denom.df,
              statname=statistic, method=method, varNames=vn)
  class(ans) <- "citest"
  ans
}

 
