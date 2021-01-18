#' Calculate the Log-Rank-Test very fast
#' @param groupa vector of group a's survival times
#' @param groupb vector of group b's survival times
#' @param groupacensored vector of censored information of group a's survival times
#' @param groupbcensored vector of censored information of group b's survival times
#' @param onlyz (optional) calculate only z-statistic
#' @return chi2 statistic, z-statistic, p-value
#' @examples
#' T1 <- c(6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25, 32, 32, 34, 35)
#' E1 <- c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
#' T2 <- c(1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23)
#' E2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' logrank_test(T1, T2, E1, E2)
#' #1.679294e+01 -4.097919e+00, 4.168809e-05
#' @export
#' @import Rcpp
#' @useDynLib fastlogranktest, .registration = TRUE
logrank_test <- function(groupa, groupb, groupacensored, groupbcensored, onlyz=FALSE){
    logrank_instance(groupa, groupb, groupacensored, groupbcensored, onlyz)
}
#' Calculate multiple Log-Rank-Tests very fast
#' @param groupas list of vectors of groupa's survival times
#' @param groupbs list of vectors of groupb's survival times
#' @param groupacensoreds list of vectors of censored information of groupa's survival times
#' @param groupbcensoreds list of vectors of censored information of groupb's survival times
#' @param threadnumber (optional) set the number of threads used for this function
#' @param onlyz (optional) calculate only z-statistic
#' @return vector of chi2 statistic, z-statistic, p-value (same order as input)
#' @examples
#' T1 <- c(6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25, 32, 32, 34, 35)
#' E1 <- c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
#' T2 <- c(1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23)
#' E2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' t1s<-list(T1, T1, T1)
#' e1s<-list(E1, E1, E1)
#' t2s<-list(T2, T2, T2)
#' e2s<-list(E2, E2, E2)
#' multi_logrank_test(t1s, t2s, e1s, e2s)
#' #1.679294e+01 -4.097919e+00 4.168809e-05 1.679294e+01 -4.097919e+00 4.168809e-05
#' #1.679294e+01 -4.097919e+00 4.168809e-05
#' @export
#' @import Rcpp
#' @useDynLib fastlogranktest, .registration = TRUE
multi_logrank_test <- function(groupas, groupbs, groupacensoreds, groupbcensoreds, threadnumber=NULL, onlyz=FALSE){
  if(is.null(threadnumber)){
      cpu_parallel_logrank(groupas, groupbs, groupacensoreds, groupbcensoreds, onlyz)
  }
  else{
      cpu_parallel_logrank1(groupas, groupbs, groupacensoreds, groupbcensoreds, threadnumber, onlyz)
  }
}
