################################################################################
### File: pseudoranks.R
### Description: Function to calculate mid-pseudo-ranks
###
################################################################################


#' Calculation of Pseudo-Ranks
#'
#' @description Calculation of (mid) pseudo-ranks of a sample. In case of ties (i.e. equal values), the average of min pseudo-rank and max-pseudor-rank are taken (similar to rank with ties.method="average").
#' @param data numerical vector
#' @param group vector coding for the groups
#' @return Returns a numerical vector containing the pseudo-ranks
#' @keywords internal
recursiveCalculation <- function(data, group, na.last, ties.method) {

  stopifnot(is.numeric(data), is.factor(group))
  n <- table(group)
  

 # case: missing values in the data
 if(sum(is.na(data)) > 0) {
   nas <- which(is.na(data))
    
   # variant 1: remove NAs
   if(is.na(na.last)) {
     data <- data[-nas]
     group <- droplevels(group[-nas])
     # group sizes need to be adjusted because NAs were removed
     n <- table(group)
   }
   # variant 2: keep NAs and put them last
   else if(na.last == TRUE) {
     m <- max(data, na.rm = TRUE)
     data[nas] <- (m+1):(m+length(nas))
   }
   # variant 3: keep NAs and put them first
   else if(na.last == FALSE) {
     m <- min(data, na.rm = TRUE)
     data[nas] <- (m-1):(m-length(nas))
   } 
  }
    
  ord <- .Call(`_pseudorank_order_vec`, data) + 1
  data_sorted <- data[ord]
  sortback <- match(data, data_sorted)
  if(ties.method == "average") {
    return(.Call(`_pseudorank_psrankCpp`, data_sorted, group[ord], n)[sortback])
  } else if(ties.method == "min") {
    return(.Call(`_pseudorank_psrankMinCpp`, data_sorted, group[ord], n)[sortback])
  } else {
    return(.Call(`_pseudorank_psrankMaxCpp`, data_sorted, group[ord], n)[sortback])
  }
}