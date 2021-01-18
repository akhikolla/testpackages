

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for non-parametric analysis on a biv.rec object for joint cdf
#'
#' @description
#' This function calculates the joint CDF for bivariate recurrent events. Called from bivrecNP(). No user interface.
#' @param fit_data Passed from bivrecNP().
#' @param u Passed from bivrecNP().
#' @param ai Passed from bivrecNP().
#' @param CI Passed from bivrecNP().
#'
#' @return A dataframe with the joint CDF
#'
#' @noRd
#' @keywords internal
#' @useDynLib BivRec

nonparam_cdf <- function(fit_data, u, ai, CI) {

  n <- fit_data$n
  m <- fit_data$m
  mc <- fit_data$mc
  nd <- fit_data$nd
  tot <- fit_data$tot
  gap <- fit_data$gap
  event <- fit_data$event
  markvar1 <- fit_data$markvar1
  markvar2 <- fit_data$markvar2
  udt <- fit_data$udt
  ctime <- fit_data$ctime
  ucen <- fit_data$ucen
  r <- fit_data$r
  d <- fit_data$d
  sest <- fit_data$sest
  Fest <- fit_data$Fest
  var <- fit_data$var
  prob <- fit_data$prob
  std <- fit_data$std
  gtime <- fit_data$gtime
  cen <- fit_data$cen
  mark1 <- fit_data$mark1
  mark2 <- fit_data$mark2

  estcdf <- list()

  for (u_count in 1:nrow(u)) {
    u1 <- u[u_count, 1]
    u2 <- u[u_count, 2]


    tmpindex <-sum(as.integer(udt<=(u1+u2)))  ### index ORINALLY PART OF BIVGAP FUNCTION
    if (tmpindex==0) {
      temp <- data.frame(u1, u2, prob=0, std=0)
      rownames(temp) <- "1"
      estcdf[[u_count]] <- temp
    } else {
      estimates <- r_bivrecur(n, gtime, ctime, mc, m,
                              cen, ucen, nd, udt, tot, gap, event,
                              r, d, sest, var, markvar1, markvar2,
                              mark1, mark2, u1, u2, Fest, tmpindex, prob, std)
      estcdf[[u_count]] <- data.frame(u1, u2, prob=estimates[1], std=estimates[2])
    }
  }

  out1 <- data.frame(matrix(unlist(estcdf), nrow=nrow(u), byrow=T))

  conf_lev = 1 - ((1-CI)/2)
  out1$lower <- out1[,3] - qnorm(conf_lev)*out1[,4]
  out1$upper <- out1[,3] + qnorm(conf_lev)*out1[,4]
  out1$lower[which(out1$lower<0)] <- 0
  out1$upper[which(out1$upper>1)] <- 1

  lowstring <-   lowstring <- paste("Lower", substr(as.character(CI), 2,4), sep=" ")
  upstring <- paste("Upper", substr(as.character(CI), 2,4), sep=" ")
  colnames(out1) <- c("x", "y", "Joint Probability", "SE", lowstring, upstring)

  return(cdf=out1)

}
