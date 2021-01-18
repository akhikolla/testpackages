#' Imputation by Simple Rules
#'
#' One of the most simplest ways to fill in the missing entries is
#' to apply any simple rule for each variable. In this example, we provide
#' 3 options, \code{"mean"}, \code{"median"}, and \code{"random"}. It assumes
#' that every column has at least one non-missing entries in that for each column,
#' the rule is applied from the subset of non-missing values.
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param method simple rule to fill in the missing entries in a columnwise manner.
#'
#' @return a named list containing \describe{
#' \item{X}{an \eqn{(n\times p)} matrix after completion.}
#' }
#'
#' @examples
#' ## load image data of 'lena128'
#' data(lena128)
#'
#' ## transform 5% of entries into missing
#' A <- aux.rndmissing(lena128, x=0.05)
#'
#' ## apply all three methods#'
#' fill1 <- fill.simple(A, method="mean")
#' fill2 <- fill.simple(A, method="median")
#' fill3 <- fill.simple(A, method="random")
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="original")
#' image(fill1$X, col=gray((0:100)/100), axes=FALSE, main="method:mean")
#' image(fill2$X, col=gray((0:100)/100), axes=FALSE, main="method:median")
#' image(fill3$X, col=gray((0:100)/100), axes=FALSE, main="method:random")
#' par(opar)
#'
#' @references
#' \insertRef{gelman_data_2007}{filling}
#'
#' @rdname fill_simple
#' @export
fill.simple <- function(A, method=c("mean","median","random")){
  #-----------------------------------------------------------------
  ## PREPROCESSING
  #   1. data check and dimension
  X = check_data(A)
  n = nrow(X)
  p = ncol(X)
  #   2. check method
  if (missing(method)){
    method = "mean"
  } else {
    method = match.arg(method)
  }
  #   3. No Columns full of NAs are allowed
  if (check_bycol(X)==FALSE){
    stop("* fill.simple : any column full of missing entries can't be handled.")
  }

  #-----------------------------------------------------------------
  ## MAIN COMPUTATION
  for (i in 1:p){
    if (method=="mean"){
      X[,i] = simple_mean(X[,i])
    } else if (method=="median"){
      X[,i] = simple_median(X[,i])
    } else if (method=="random"){
      X[,i] = simple_random(X[,i])
    }
  }

  #-----------------------------------------------------------------
  ## RETURN OUTPUT
  result = list()
  result$X = X
  return(result)
}




#' @keywords internal
#' @noRd
simple_random <- function(a){
  missing <- is.na(a)
  n.missing <- sum(missing)
  a.obs <- a[!missing]
  imputed <- a
  imputed[missing] <- sample(a.obs,n.missing,replace=TRUE)
  return (imputed)
}
#' @keywords internal
#' @noRd
simple_mean <- function(a){
  missing <- is.na(a)
  n.missing <- sum(missing)
  a.obs <- a[!missing]
  imputed <- a
  imputed[missing] <- mean(a.obs)
  return (imputed)
}
#' @keywords internal
#' @noRd
simple_median <- function(a){
  missing <- is.na(a)
  n.missing <- sum(missing)
  a.obs <- a[!missing]
  imputed <- a
  imputed[missing] <- median(a.obs)
  return (imputed)
}
