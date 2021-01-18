###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for non-parametric analysis on a biv.rec object
#'
#' @description
#' This function calculates the marginal survival for bivariate recurrent events. Called from biv.rec.np(). No user interface.
#' @param fit_data An object that has been reformatted using the biv.rec.reformat() function. Passed from biv.rec.np().
#' @param CI Passed from biv.rec.np().
#'
#' @return A data frame with marginal survival
#'
#' @noRd
#' @keywords internal
#'

nonparam_marginal <- function(fit_data, CI) {
  n <- fit_data$n
  m <- fit_data$m
  mc <- fit_data$mc
  nd <- fit_data$nm1
  tot <- fit_data$tot
  gap <- fit_data$markvar1
  event <- fit_data$event
  udt <- fit_data$umark1
  ctime <- fit_data$ctime
  ucen <- fit_data$ucen
  gtime <- fit_data$mark1
  cen <- fit_data$cen
  r = d = sest = std = rep(0, nd)

  surv <- r_onesamp(n,gtime,ctime,mc,m,
                    cen,ucen,nd,udt,tot,gap,event,
                    r,d,sest,std)

  surv <- r_onesamp(n,gtime,ctime,mc,m,
                    cen,ucen,nd,udt,tot,gap,event,
                    r,d,sest,std)

  conf_lev = 1 - ((1-CI)/2)
  surv$lower <- surv[,2] - qnorm(conf_lev)*surv[,3]
  surv$upper <- surv[,2] + qnorm(conf_lev)*surv[,3]
  surv$lower[which(surv$lower<0)] <- 0
  surv$upper[which(surv$upper>1)] <- 1

  lowstring <- paste("Lower", substr(as.character(CI), 2,4), sep=" ")
  upstring <- paste("Upper", substr(as.character(CI), 2,4), sep=" ")
  colnames(surv) <- c("Time", "Marginal Survival", "SE", lowstring, upstring)

  return(marg_survival = surv)

}
