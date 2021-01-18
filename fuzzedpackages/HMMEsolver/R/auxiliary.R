#' @keywords internal
#' @noRd
check_matrix <- function(X,n=nrow(X),p=ncol(X)){
  cond1 = is.matrix(X)
  cond2 = all(!is.na(X))
  cond3 = all(!is.infinite(X))
  cond4 = ((n==nrow(X))&&(p==ncol(X)))
  if (cond1&&cond2&&cond3&&cond4){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
check_vector <- function(x,n=length(x)){
  cond1 = is.vector(x)
  cond2 = all(!is.na(x))
  cond3 = all(!is.infinite(x))
  cond4 = (length(x)==n)
  if (cond1&&cond2&&cond3&&cond4){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
