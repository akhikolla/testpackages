## CHECKER For Various Cases
#   0. islist_multiple     : list of multiple length
#   1. islist_sphere       : must be vector, norm 1
#   2. islist_spd          : symmetric, positive definite matrices
#   3. islist_euclidean    : only check vector or matrix
#   4. islist_stiefel      : stiefel and grassman both can be used
#
#   6. islist_rotation     : rotation group SO(n)

# 0. islist_multiple ------------------------------------------------------
#' @keywords internal
#' @noRd
islist_multiple <- function(data){
  cond1 = is.list(data)
  cond2 = (length(data)>=2)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# 1. islist_sphere --------------------------------------------------------
#' @keywords internal
#' @noRd
islist_sphere <- function(data){
  single_checker <- function(x){
    y = as.vector(x)
    if (abs(sum(y^2)-1)<sqrt(.Machine$double.eps)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  
  if (is.null(nrow(data[[1]]))){
    cond1 = (all(unlist(lapply(data, is.vector))==TRUE)) # must be all vectors
    cond2 = all(unlist(lapply(data, single_checker))==TRUE)
    cond3 = (length(unique(unlist(lapply(data, length))))==1)
    if (cond1&&cond2&&cond3){
      return(TRUE)
    } else {
      return(FALSE)
    }  
  } else {
    return(FALSE)
  }
}

# 2. islist_spd -----------------------------------------------------------
#' @keywords internal
#' @noRd
islist_spd <- function(data){
  if (!is.null(nrow(data[[1]]))){
    
    cond1.diag = (all(unlist(lapply(data, is.matrix))==TRUE))                           # must be all matrices
    if (!cond1.diag){
      stop("* riemfactory : 'spd' requires all data to be square matrices.")
    }
    cond1 = cond1.diag
    cond2.diag = (all(unlist(lapply(data, isSymmetric.matrix))==TRUE))                  # all symmetric matrices
    if (!cond2.diag){
      stop("* riemfactory : 'spd' requires all matrices be symmetric.")
    }
    cond2 = cond2.diag
    cond3 = (length(unique(unlist(lapply(data, nrow))))==1)                        # all have same sizes
    
    vec4  = rep(0,length(data))
    for (i in 1:length(data)){
      vec4[i] = min(eigen(data[[i]])$values)
    }
    cond4 = (all(vec4>=0))
    if (!cond4){
      stop("* riemfactory : 'spd' takes positive semi-definite matrices only.")
    }
    if (cond1&&cond2&&cond3&&cond4){
      return(TRUE)
    } else {
      return(FALSE)
    }  
  } else {
    return(FALSE)
  }
}

# 3. islist_euclidean -----------------------------------------------------
#' @keywords internal
#' @noRd
islist_euclidean <- function(data){
  if (is.vector(data[[1]])){
    cond1 = (all(unlist(lapply(data, is.vector)))==TRUE)
    cond2 = (length(unique(unlist(lapply(data, length))))==1)
    if (cond1&&cond2){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else if (is.matrix(data[[1]])){
    cond1 = (length(unique(unlist(lapply(data, nrow))))==1)
    cond2 = (length(unique(unlist(lapply(data, ncol))))==1)
    cond3 = (all(unlist(lapply(data, is.matrix)))==TRUE)
    if (cond1&&cond2&&cond3){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else{
    return(FALSE)
  }
}


# 4. islist_stiefel -------------------------------------------------------
#' @keywords internal
#' @noRd
islist_stiefel <- function(data){
  if (!is.null(nrow(data[[1]]))){ # if a matrix
    cond1 = (length(unique(unlist(lapply(data, nrow))))==1)
    if (!cond1){
      stop("* riemfactory : 'stiefel' requires all data to have same number of rows.")
    }
    cond2 = (length(unique(unlist(lapply(data, ncol))))==1)
    if (!cond2){
      stop("* riemfactory : 'stiefel' requires all data to have same number of columns.")
    }
    cond3 = all(unlist(lapply(data, function(X){(norm(t(X)%*%X - diag(ncol(X))) < sqrt(length(data)*.Machine$double.eps))}))==TRUE)
    if (!cond3){
      stop("* riemfacotyr : 'sitefel' matrix inner product is not equal to identity matrix.")
    }
    if (cond1&&cond2&&cond3){
      return(TRUE)
    } else {
      return(FALSE)
    }  
  } else {
    return(FALSE)
  }
}


# 6. islist_rotation     : rotation group SO(n) ---------------------------
#' @keywords internal
#' @noRd
islist_rotation <- function(data){
  if (!is.null(nrow(data[[1]]))){ # if a matrix
    cond1 = (length(unique(unlist(lapply(data, nrow))))==1)
    cond2 = (length(unique(unlist(lapply(data, ncol))))==1)
    cond3 = (nrow(data[[1]])==ncol(data[[2]]))
    n     = nrow(data[[1]])
    thr   = sqrt(123*.Machine$double.eps)
    cond4 = all(unlist(lapply(data, function(X){norm(diag(n)-(t(X)%*%X),"f")<thr}))==TRUE) # R^T R  = I
    cond5 = all(unlist(lapply(data, function(X){abs(det(X)-1) < thr}))==TRUE)              # det(R) = 1
    if (cond1&&cond2&&cond3&&cond4&&cond5){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}


