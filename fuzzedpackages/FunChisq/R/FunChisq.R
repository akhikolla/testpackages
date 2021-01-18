# FunChisq.R -- statistical tests for model-free functional dependency
#
# YZ, HZ, MS
# Modified: Feb 8, 2015; June 25, 2015
# Jan 22, 2016 MS:
#    1. introduced the type argument to indicate whether alternative hypothesis is functional
#       or non-constant function
#    2. introduced the log.p argument to obtain the log of p-value
# Jan 23, 2016 MS:
#    Added a return of an estimate of functional index between 0 and 1. It is asymetrical and
#      different from Cramer's V.
# Jan 30, 2016 MS:
#    Added a new argument "index.kind" to specify the function index kind: "unconditional"
#      or "conditional" given marginal of Y
# Feb 5, 2016 MS:
#    Renamed the "type" argument to "alternative"
# Feb 7, 2016 MS:
#   Handled the special case when the degrees of freedom for the normalized FunChisq are zero
# May 3, 2016 MS:
#   Listed all arguments
# May 29, 2017 HZ:
#   Added "simulate.p.value" method to calculate p-value with simulated distribution
# Feb 08, 2018 HZ:
#   Added exact functional test option "exact.mode.bound"
# Oct 24, 2018 MS:
#   Default to (conditional) function index. The previous default was unconditional function
#   index.

fun.chisq.test <- function (
  x,
  method=c("fchisq", "nfchisq",
           "exact", "exact.qp", "exact.dp", "exact.dqp",
           "default","normalized", "simulate.p.value"),
  alternative=c("non-constant", "all"), log.p=FALSE,
  index.kind=c("conditional", "unconditional"
#               , "fix.row.sums", "fix.column.sums",
#               "fix.marginal.sums"
               ),
  simulate.nruns = 2000,
  exact.mode.bound = TRUE
)
{
  if(!is.matrix(x) && !is.data.frame(x)) stop("input x must be matrix or data frame\n")

  method <- match.arg(method)
  if(method == "default") {
    warning(paste0("method=\"", method, "\" is deprecated. Use \"fchisq\" instead."))
    method <- "fchisq"
  } else if(method == "normalized") {
    warning(paste0("method=\"", method, "\" is deprecated. Use \"nfchisq\" instead."))
    method <- "nfchisq"
  }

  alternative <- match.arg(alternative)
  index.kind <- match.arg(index.kind)

  row.sums <- rowSums(x)
  col.sums <- colSums(x)

  n <- sum(col.sums)

  if(0) {
    row.chisq.sum <- sum(apply(x, 1,
                               function(v){
                                 row.sum <- sum(v)
                                 expected <- row.sum / length(v)
                                 if(row.sum>0) sum( (v - expected)^2 / expected)
                                 else 0
                               }))
  } else { # MS June 13, 2017
    row.chisq.sum <- sum(sapply(seq(nrow(x)),
                                function(i){
                                  if(row.sums[i]>0) sum( x[i, ]^2 ) / row.sums[i]
                                  else 0
                                })) * ncol(x) - n
  }

  if(alternative == "non-constant") { # non-constant functional chi-squared
    # col.expected <- n / ncol(x)
    # col.chisq <- sum((col.sums - col.expected)^2 / col.expected)

    col.chisq <- sum(col.sums^2) / n * ncol(x) - n

    fun.chisq <- row.chisq.sum - col.chisq
    df <- nrow(x) * (ncol(x) - 1) - (ncol(x) - 1)

    # Version 1: Current. Added June 11, 2017. Bound is maximum reachable
    if(index.kind == "unconditional") { # unconditional functional index
      max.fun.chisq <- n * ncol(x) * (1 - 1.0 / min(dim(x)))
      estimate.label <- "unconditional function index xi.f"
    } else if(index.kind == "conditional" || index.kind == "fix.column.sums") {
      # functional index given the column marginal

      # max.fun.chisq <- sum(sort(col.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
      #  ncol(x) - n - col.chisq

      max.fun.chisq <- n * ncol(x) - n - col.chisq

      estimate.label <- "function index xi.f" # "function index xi.f conditioned on column sums"

    } else if(index.kind == "fix.row.sums") { # Y-conditional functional index

      max.fun.chisq <- sum(sort(row.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
        ncol(x) - n - col.chisq
      estimate.label <- "function index xi.f conditioned on row sums (inaccurate)"

    } else if(index.kind == "fix.marginal.sums") { # XY-conditional functional index

      max.fc.given.col.sums <- sum(sort(col.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
        ncol(x) - n - col.chisq

      max.fc.given.row.sums <- sum(sort(row.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
        ncol(x) - n - col.chisq

      max.fun.chisq <- min(c(max.fc.given.col.sums, max.fc.given.row.sums))

      estimate.label <- "function index xi.f conditioned on row and column sums (inaccurate)"
    }

    # Version 0: bound is not reachable when nrow(x) < ncol(x)
    else if(index.kind == "conditional-version-0") { # functional index given the column marginal
      max.fun.chisq <- n * (ncol(x) - 1) - col.chisq
      estimate.label <- "conditional function index xi.f"
    } else if(index.kind == "unconditional-version-0") { # unconditional functional index
      max.fun.chisq <- n * (ncol(x) - 1)
      estimate.label <- "function index xi.f"
    } else {
      stop("unrecognized function index kind ", index.kind)
    }

  } else if(alternative == "all") { # functional chi-squared

    fun.chisq <- row.chisq.sum
    df <- nrow(x) * (ncol(x) - 1)

    if(index.kind == "unconditional") {
      max.fun.chisq <- n * (ncol(x) - 1)
      estimate.label <- "function index xi.f"

    } else if(index.kind == "conditional" || index.kind == "fix.column.sums") {
      # MS. Added June 11, 2017
      # max.fun.chisq <- sum(sort(col.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
      #  ncol(x) - n

      max.fun.chisq <- n * ncol(x) - n

      estimate.label <- "function index xi.f conditioned on column sums"

    } else if(index.kind == "fix.row.sums") { # Y-conditional functional index

      max.fun.chisq <- sum(sort(row.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
        ncol(x) - n
      estimate.label <- "function index xi.f conditioned on row sums"

    } else if(index.kind == "fix.marginal.sums") { # XY-conditional functional index

      max.fc.given.col.sums <- sum(sort(col.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
        ncol(x) - n

      max.fc.given.row.sums <- sum(sort(row.sums, decreasing=TRUE)[seq(min(dim(x)))]) *
        ncol(x) - n

      max.fun.chisq <- min(c(max.fc.given.col.sums, max.fc.given.row.sums))
      estimate.label <- "function index xi.f conditioned on row and column sums"
    }

  } else {
    stop("unknown alternative ", alternative)
  }

  estimate.label <- paste(alternative, estimate.label)

  if(max.fun.chisq > 0) {
    estimate <- sqrt(abs(fun.chisq) / max.fun.chisq)
  } else {
    estimate <- 0
  }

  names(estimate) <- estimate.label

  DNAME <- deparse(substitute(x))

  if(method=="fchisq") {
    method.text <- "Functional chi-squared test"
    names(fun.chisq) <- "statistic"
    names(df) <- "parameter"
    p.value <- pchisq( fun.chisq, df = df, lower.tail=FALSE, log.p=log.p )
    return(structure(list( statistic=fun.chisq, parameter=df, p.value=p.value,
                           estimate = estimate, data.name=DNAME,
                           method = method.text),
                     class = "htest"))

  } else if(method=="nfchisq") {
    method.text <- "Normalized functional chi-squared test"
    if(df > 0) {
      normalized <- as.numeric((fun.chisq-df)/sqrt(2*df))
    } else {
      normalized <- -Inf
    }
    p.value <- pnorm( normalized, lower.tail=FALSE, log.p=log.p )
    names(normalized) <- "statistic"
    names(df) <- "parameter"
    return(structure(list(statistic = normalized, parameter = df, p.value = p.value,
                          estimate = estimate, data.name = DNAME, method = method.text),
                     class = "htest"))
  } else if(method=="simulate.p.value"){
    method.text <- "Functional chi-squared test with simulated p value"
    names(fun.chisq) <- "statistic"
    names(df) <- "parameter"
    #p.value <- pchisq( fun.chisq, df = df, lower.tail=FALSE, log.p=log.p )
    p.value <- simulate.p.value(x, simulate.nruns)

    return(structure(list( statistic=fun.chisq, parameter=df, p.value=p.value,
                           estimate = estimate, data.name=DNAME,
                           method = method.text),
                     class = "htest"))
  } else if(method=="exact.qp") {

    method.text <- "Exact functional test"
    if(sum(x%%1!=0)>=1) { # Check whether numbers in x are all integers
      stop("ERROR: Exact test requires integer contingency tables!", call. = TRUE)
    }
    ####

    if((sum(x) <= 200 || sum(x)/nrow(x)/ncol(x) <=5)
       && nrow(x)<=10 && ncol(x)<=10) {
      # p.value <- exact.functional.test(x)
      p.value <- ExactFunctionalTest(x, exact.mode.bound)
      if(log.p) p.value <- log(p.value)
      names(fun.chisq) <- "statistic"
      return(structure(list(statistic = fun.chisq, p.value = p.value, estimate = estimate,
                            data.name = DNAME, method = method.text),
                       class = "htest"))
    } else {
      warning("Asymptotic test is used in place of exact test to save time.\n")
      return(fun.chisq.test(x, method="fchisq", alternative=alternative, log.p=log.p,
                            index.kind=index.kind))
    }
  } else if(method=="exact.dp") {

    method.text <- "Exact functional test"
    if(sum(x%%1!=0)>=1) { # Check whether numbers in x are all integers
      stop("ERROR: Exact test requires integer contingency tables!", call. = TRUE)
    }
    ####

    if((sum(x) <= 200 || sum(x)/nrow(x)/ncol(x) <=5)
       && nrow(x)<=10 && ncol(x)<=10) {

      p.value <- EFTDP(x)
      if(log.p) p.value <- log(p.value)
      names(fun.chisq) <- "statistic"
      return(structure(list(statistic = fun.chisq, p.value = p.value, estimate = estimate,
                            data.name = DNAME, method = method.text),
                       class = "htest"))
    } else {
      return(fun.chisq.test(x, method="fchisq", alternative=alternative, log.p=log.p,
                            index.kind=index.kind))
    }
    ####
  } else if(method=="exact" || method=="exact.dqp") {

    method.text <- "Exact functional test"
    if(sum(x%%1!=0)>=1) { # Check whether numbers in x are all integers
      stop("ERROR: Exact test requires integer contingency tables!", call. = TRUE)
    }
    ####

    # if((sum(x) <= 200 || sum(x)/nrow(x)/ncol(x) <=5)
    #   && nrow(x)<=10 && ncol(x)<=10)
    {

      p.value <- EFTDQP(x)
      if(log.p) p.value <- log(p.value)
      names(fun.chisq) <- "statistic"
      return(structure(list(statistic = fun.chisq, p.value = p.value, estimate = estimate,
                            data.name = DNAME, method = method.text),
                       class = "htest"))
    } ## else {
      ## return(fun.chisq.test(x, method="fchisq", alternative=alternative, log.p=log.p,
      ##                      index.kind=index.kind))
    ## }
    ####
  }
  else {
    stop("ERROR: unrecognized method argument", method)
  }
}
