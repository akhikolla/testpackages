# AUXILIARY FUNCTIONS -----------------------------------------------------
# 01. check_data  : common.
# 02. check_bycol : maybe only for simple methods
# 03. check_na    : it is inclusive of "is.na" and "is.nan"

# 01. check_data ----------------------------------------------------------
# NA is okay, but Inf case, send a warning
#' @keywords internal
#' @noRd
check_data <- function(X){
  # 1. should NOT be a constant
  if (length(as.vector(X))==1){
    stop("* filling : input data should not be a single number.")
  }
  # 2. vector into matrix (1-by-length(X))
  if (is.vector(X)){
    X = matrix(X,ncol=length(X))
  }
  # 3. Inf case : send warning
  if (any(is.infinite(X))){
    stop("* filling : Inf or -Inf values are not allowed. Missing entries should be NAs.")
  }
  return(X)
}

# 02. check_bycol ---------------------------------------------------------
# No column is permitted where every entry is NA
# Returns FALSE if it is NOT CLEAN, i.e., there exists one column at least of all NAs
#         TRUE  if it is CLEAN
#' @keywords internal
#' @noRd
check_bycol <- function(X){
  fun <- function(x){(sum(is.na(x))==length(x))||(sum(is.nan(x))==length(x))}
  res <- apply(X, 2, fun)
  if ((sum(res))>0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}


# 03. check_na ------------------------------------------------------------
#' @keywords internal
#' @noRd
check_na <- function(X){
  if ((is.na(X))||(is.nan(X))){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
