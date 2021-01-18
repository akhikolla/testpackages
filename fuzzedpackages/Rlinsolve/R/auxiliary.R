
# Auxiliary Functions -----------------------------------------------------
# Aux 1:: aux.is.dd(w/ doc) : check if diagonally dominant
# Aux 2:: aux.is.psd    : check positive semidefinite
# Aux 3:: aux.is.sparse : sparseformats = c("dgCMatrix","dtCMatrix","dsCMatrix")



#   -----------------------------------------------------------------------
# Aux 3:: aux.is.sparse : check whether one of the following
#' @keywords internal
#' @noRd
aux.is.sparse <- function(AA){
  sparseformats = c("dgCMatrix","dtCMatrix","dsCMatrix")
  for (i in 1:3){
    if (inherits(AA, sparseformats[i])){
      return(TRUE)
    }
  }
  return(FALSE)
}

#  ------------------------------------------------------------------------
# Aux 1:: aux.is.dd
#   "sdd" : strictly
#   "wdd" : weakly
#   FALSE : not at all
#' @keywords internal
#' @noRd
aux.is.dd <- function(A){
  absA = abs(A)
  # 1-1. separate terms
  diagA = 2*(diag(absA))
  offdA = colSums(absA)
  # 1-2. logic
  if (all(diagA>offdA)){
    res = "sdd"
  } else if ((all(diagA>=offdA))&&(any(diagA==offdA))){
    res = "wdd"
  } else {
    res = FALSE
  }
  return(res)
}

#  ------------------------------------------------------------------------
# Aux 2:: aux.is.psd
#' Positive Semidefiniteness
#' PD, PSD, or FALSE
#'
#' @keywords internal
#' @noRd
aux.is.psd <- function(A){
  # get eigenvalues
  eigs = eigen(A, only.values = TRUE)

  # PD, PSD, or FALSE
  if (all(eigs>0)){res = "PD"}
  else if ((all(eigs>=0))&&(any(eigs>0))){res = "PSD"}
  else {res = FALSE}

  # finalize
  return(res)
}
