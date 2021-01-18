np_dat <- function(dat, ai) {

  id <- dat$id
  uid <- unique(id)   # vector of unique id's
  n.uid <- length(uid)   # scalar : number of unique IDs
  event <- dat$d2 #event indicator : must always be 0 for the last obs per ID and 1 otherwise
  markvar1 <- dat$vij
  markvar2 <- dat$wij
  gap <- markvar1 + markvar2

  m.uid <- as.integer(table(id))   # vector: number of observed pairs per id/subject (m)
  max.m <- max(m.uid, na.rm=T) # scalar : maximum number of repeated observations

  ifelse (ai == 1, weight <- rep(1, n.uid), weight <- dat$ci[which(dat$epi == 1)]) #Set weights

  tot <- length(gap) # total number of observations

  ugap <- sort(unique(gap[event == 1]))   # sorted unique uncensored X_0 gap times (support points for sum)
  n.ugap <- length(ugap)   # number of unique X_0 gap times (or support points for sum)
  umark1 <- sort(unique(markvar1[event == 1]))   # sorted unique uncensored V_0 times (support points for marginal)
  n.umark1 <- length(umark1) # number of unique V_0 gap times (or support points for marginal)


  # Space holders
  r <- sest <- Fest <- rep(0, n.ugap)
  d <- matrix(0, nrow = n.ugap, ncol = 2)
  prob <- var <- std <- 0
  gtime <- cen <- mark1 <- mark2 <- matrix(0, nrow = n.uid, ncol = max.m)

  out <- list(n = n.uid, m = m.uid, mc = max.m, nd = n.ugap, tot=tot,
              gap =gap, event = event, markvar1 = markvar1, markvar2 =markvar2,
              udt = ugap,  ctime = weight, ucen = m.uid-1,
              r = r, d=d, sest = sest, Fest = Fest, var = var,
              prob = prob, std = std, gtime = gtime, cen = cen,
              mark1 = mark1, mark2 = mark2, umark1=umark1, nm1 = n.umark1)

  return(out)
}

np_fit4conditional <- function(data, ai, u1, u2){

  ### PULL INFORMATION FROM PARAMETERS TO SEND TO REFORMAT
  # identifier=xij=yij=c_indicatorY=c_indicatorX=episode=covariates=NULL
  # method <- "Non-Parametric"
  # condgx <- TRUE

  ###Send to biv.rec.reformat and complete analysis
  #new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai, condgx, data)

  my_data = na.omit(data)
  forcdf <- np_dat(dat=my_data, ai=ai)
  new_data <- list(forcdf=forcdf, refdata = my_data)
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)
  res1 <- nonparam_cdf(fit_data=new_data$forcdf, u, ai, CI=0.95)[,1:4]

  return(res1)
}

bstp <- function(seedi, ps1, ps2, x.grid, y.grid, n, refdata, ai, mintime) {
  set.seed(seedi)
  samp.id <- sample(1:n, n, replace = TRUE)
  bootdat <- NULL
  for (j in 1:n){
    temp <- refdata[which(refdata$id==samp.id[j]), ]
    bootdat <- rbind(bootdat, cbind(id = j, temp[, -1]))
  }

  joint2 <- np_fit4conditional(data=bootdat, #got rid of formula
                               ai=ai, u1=x.grid$Time[2], u2=y.grid)

  if (x.grid$Time[1] == mintime) {
    conditional <- joint2[,3] / (1-ps2)
  } else {
    joint1 <- np_fit4conditional(data=bootdat,
                                 ai=ai, u1=x.grid$Time[1], u2=y.grid) #got rid of formula
    conditional <- (joint2[,3] - joint1[,3])/(ps1 - ps2)
  }

  return(conditional)
}

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for additional non-parametric analysis of bivariate recurrent event data.
#'
#' @description
#' This function calculates the conditional cdf after estimation of joint cdf and marginal survival.  Called from Passed from bivrecNP(). No user interface.
#' @param res List with joint.cdf and marginal.survival. Passed from bivrecNP()
#' @param given.interval is a 1x2 vector indicating an interval for the first gap time to estimate the cdf of second gap time. Passed from Passed from bivrecNP()
#' @param CI confidence level. Passed from bivrecNP()
#' @param condiplot a logical value. Passed from bivrecNP()
#' @param yij Passed from bivrecNP()
#'
#' @return A data frame with the conditional CDF for the given an interval of the first gap time and corresponding plot.
#' @importFrom stats sd
#'
#' @noRd
#' @keywords internal
#'

nonparam_conditional <- function(res, given.interval, CI, yij) {

  ####Extract items from results
  marginal <- res$marginal_survival #marg result (res2)
  marginal$rounded <- round(marginal[,2], digits=2)
  ai <- res$ai
  new_data <- res$new_data #this is essentially the "fit_data"
  margdata <- new_data$formarg
  refdata <- new_data$refdata
  n <- margdata$n

  x.grid <- marginal[which(marginal$Time<=given.interval[2]), ] #this uses marginal result and is a param for bstp
  x.grid <- x.grid[which(x.grid$Time>=given.interval[1]), ]
  x.grid$diffs <- x.grid$rounded - x.grid$rounded[nrow(x.grid)]
  if (length(which(x.grid$diffs >= 0.1))==0) {
    print("Cannot estimate conditional cdf, given.interval is too narrow")
    stop()
  }
  x.grid <- x.grid[c(1, max(which(x.grid$diffs >= 0.05)), nrow(x.grid)),]
  y.grid <- seq(min(yij), max(yij), length.out = 200)
  ps1 <- x.grid[1,2]
  ps2 <- x.grid[3,2]

  B = ifelse(CI==0.99, 200, 100)
  cond.prob <- matrix(rep(NA, length(y.grid)*B), ncol=B)
  colnames(cond.prob) = seq(1,B,1)
  print(paste("Estimating conditional cdf with ", CI*100, "% confidence interval using ", B, " bootstrap samples", sep=""))

  for (i in 1:B) {
    #print(paste("Sample", i, sep = " "))
    cond.prob[,i] <- bstp(seedi=i, ps1, ps2, x.grid, y.grid, n, refdata, ai, mintime = min(marginal$Time))
  }

  conf.lev = 1 - ((1-CI)/2)
  bootstrapCIs <- apply(cond.prob, 1, function(x) c(mean(x), sd(x), sort(x)[(1-conf.lev)*B], sort(x)[conf.lev*B]))
  cond <- round(data.frame(y.grid, bootstrapCIs[1,], bootstrapCIs[2,], bootstrapCIs[3,], bootstrapCIs[4,]), digits = 4)
  low.string <- paste("Bootstrap Lower", substr(as.character(CI), 2,4), sep=" ")
  up.string <- paste("Bootstrap Upper", substr(as.character(CI), 2,4), sep=" ")
  colnames(cond) <- c("Time", "Conditional Probability" , " Bootstrap SE", low.string, up.string)

  flat.ind <- which(cond[,5]>=1.001)
  if (length(flat.ind)!=0) {cond[flat.ind, 2:5] <- cond[(min(flat.ind)-1), 2:5]}
  #if (condiplot == TRUE) {
  #plot(cond$Time, cond[,5], type="l", lty = 2, xlab = "Type II Gap Times (y)", ylab = "Conditional Probability",
  #xlim=c(0, round(max(y.grid), digits=1)),
  #ylim=c(0, round(max(cond[,5]), digits=1)),
  #main=substitute(
  #paste("P(", Y^0 <= y, "|", X^0 %in% "[", gi1, ",", gi2, "])"),
  #list(gi1 = given.interval[1], gi2 = given.interval[2]))
  #)
  #graphics::lines(cond$Time, cond[,4], lty = 2)
  #graphics::lines(cond$Time, cond$Conditional.Probability,lty = 1)

  #}

  cond[, 4:5] <- round(cond[,4:5], digits = 2)
  # cond$xgrid=x.grid
  # cond$ygrid=y.grid
  # cond$data=data
  # cond$cond.prob
  return(list(conditional=cond, ygrid=y.grid))
  #return(list(conditional=cond,xgrid=x.grid,ygrid=y.grid,data=data,condprob=cond.prob,yij=yij))
}
