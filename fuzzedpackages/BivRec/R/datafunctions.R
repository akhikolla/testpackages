#                 m.dat, np.dat FUNCTIONS                                      #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson and Aparajita Sur (March, 2019)       #
# Received from Chihyun Lee (January, 2018)                                    #
#______________________________________________________________________________#

#' A function to reformat data for Lee fit
#'
#' @description
#' This function reformats data for Lee fit
#'
#' @importFrom stats na.omit
#'
#' @param dat A data frame
#'
#' @noRd
#' @keywords internal
#'
mdat <- function(dat) {
  n=length(unique(dat$id))
  mc=max(dat$epi)-1

  g1dat=cbind(dat[dat$epi==1,]$xij,1-dat[dat$epi==1,]$d1)
  g2dat=cbind(dat[dat$epi==1,]$zij,1-dat[dat$epi==1,]$d2)
  l1=max(g1dat[g1dat[,2]==0,1])-(1e-07)
  l2=max(g2dat[g2dat[,2]==0,1])-(1e-07)
  g1surv=survival::survfit(Surv(g1dat[,1],g1dat[,2])~1)
  g2surv=survival::survfit(Surv(g2dat[,1],g2dat[,2])~1)

  xmat=ymat=zmat=delta1=delta2=g1mat=g2mat=matrix(0,n,mc,byrow=TRUE)
  mstar=ctime=NULL
  for (i in 1:n) {
    tmp=dat[dat$id==i,]
    tmp.mstar=ifelse(nrow(tmp)==1,1,nrow(tmp)-1)
    mstar=c(mstar,tmp.mstar)
    ctime=c(ctime,tmp$ci[1])

    xmat[i,1:tmp.mstar]=tmp$xij[1:tmp.mstar]
    ymat[i,1:tmp.mstar]=tmp$yij[1:tmp.mstar]
    zmat[i,1:tmp.mstar]=tmp$zij[1:tmp.mstar]
    delta1[i,1:tmp.mstar]=tmp$d1[1:tmp.mstar]
    delta2[i,1:tmp.mstar]=tmp$d2[1:tmp.mstar]

    g1mat[i,1:tmp.mstar]=sapply(xmat[i,1:tmp.mstar],function(x) summary(g1surv,times=min(x,l1))$surv)
    g2mat[i,1:tmp.mstar]=sapply(zmat[i,1:tmp.mstar],function(x) summary(g2surv,times=min(x,l2))$surv)
  }

  cumh1=cumsum(g1surv$n.event/g1surv$n.risk)
  cumh2=cumsum(g2surv$n.event/g2surv$n.risk)
  l1mat=cbind(g1surv$time,diff(c(cumh1,tail(cumh1,1))),g1surv$surv)
  l2mat=cbind(g2surv$time,diff(c(cumh2,tail(cumh2,1))),g2surv$surv)

  out=list(n=n,mc=mc,xmat=xmat,ymat=ymat,zmat=zmat,delta1=delta1,delta2=delta2,g1mat=g1mat,g2mat=g2mat,l1=l1,l2=l2,l1mat=l1mat,l2mat=l2mat, mstar=mstar,ctime=ctime)
  return(out)
}

#####Reformat data set for non-parametric analysis

# dat : a data.frame including
#     1)id numbers, 2)orders of episodes, 3)first gap time, 4)second gap time
#     5)censoring times, 6) censoring indicators in each column
# ai: a non-negative function of censoring time

#' A function to reformat data for all non-parametric analysis
#'
#' @description
#' This function reformat data for all non-parametric analysis
#'
#' @importFrom stats na.omit
#'
#' @param dat a data frame including id numbers, order of episodes, first gap time, second gap time,
#' censoring times, and censoring indicators in each column
#' @param ai a non-negative function of censoring time
#'
#' @noRd
#' @keywords internal
#'
#'
np.dat <- function(dat, ai) {

  id <- dat$id
  uid <- unique(id)   # vector of unique id's
  n.uid <- length(uid)   # scalar : number of unique IDs
  event <- dat$d2 #event indicator : must always be 0 for the last obs per ID and 1 otherwise
  markvar1 <- dat$vij #gap times of type 1
  markvar2 <- dat$wij #gap times of type 2
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

#' A function to reformat data for marginal portion of nonparametric analysis
#'
#' @description
#' This function reformats data for marginal portion of nonparametric analysis
#'
#' @importFrom stats na.omit
#'
#' @param dat a data frame
#'
#' @noRd
#' @keywords internal
#'
#'
formarginal <- function(dat){

  mdata <- tmp <- NULL
  freq <-cumsum(c(0,table(dat[,1])))

  for (i in 1:(length(freq)-1)){
    tmp<- dat[(freq[i]+1):freq[i+1], -c(5,7)]
    if (nrow(tmp)==1){
      if(tmp$wij>0){
        mdata <- rbind(mdata,
                       c(id=tmp$id, vij=tmp$vij, wij=tmp$wij, d2=1, epi=tmp$epi, ci=tmp$ci),
                       c(id=tmp$id, vij=0, wij=0, d2=0, epi=2, ci=tmp$ci))
      } else{
        mdata<-rbind(mdata, c(id=tmp$id, vij=tmp$vij, wij=tmp$wij, d2=0, epi=tmp$epi, ci=tmp$ci))
      }
    } else{mdata<-rbind(mdata, tmp)}
  }
  return(as.data.frame(mdata))
}
