###########################################################################
############## FUNCTIONS FOR REFERENCE BY MAIN - NOT FOR USER #############
###########################################################################

#           m.dat.chang1, all RE, v.est1 and sd.estpar1 FUNCTIONS              #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (June, 2018)                          #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________

######################
##-----reformat dataset
m.dat.chang1=function(dat,beta) {
  n = length(unique(dat$id))
  beta1 = beta[1]
  beta2 = beta[2]
  maxb = apply(cbind(beta1,beta2), 1, max)
  amat = dat$a1
  dat$txij = dat$xij*exp(-amat*beta1)
  dat$tzij = dat$txij + dat$yij*exp(-amat*beta2)
  dat$tci = dat$ci*exp(-amat*maxb)

  all.t.xij = all.t.zij = all.t.d1 = all.t.d2 = all.mstar = all.a = NULL

  for (i in unique(dat$id)) {
    tmp=dat[dat$id==i,]

    t.xij=min(tmp$txij[1],tmp$tci[1])
    t.zij=min(tmp$tzij[1],tmp$tci[1])
    t.d1=(t.xij<tmp$tci[1])
    t.d2=(t.zij<tmp$tci[1])

    if (nrow(tmp)>1) {
      td1=t.d1
      td2=t.d2
      j=2
      while (td1==1 & td2==1 & j<=nrow(tmp)) {
        tsum=sum(t.zij[1:(j-1)])
        txij=min(tmp$txij[j],tmp$tci[j]-tsum)
        tzij=min(tmp$tzij[j],tmp$tci[j]-tsum)
        td1=(txij+tsum)<tmp$tci[j]
        td2=(tzij+tsum)<tmp$tci[j]
        if (td1==1 & td2==1) {
          t.xij=c(t.xij,txij)
          t.zij=c(t.zij,tzij)
          t.d1=c(t.d1,td1)
          t.d2=c(t.d2,td2)
        }
        j=j+1
      }
    }
    all.t.xij = c(all.t.xij,t.xij)
    all.t.zij = c(all.t.zij,t.zij)
    all.t.d1 = c(all.t.d1,t.d1)
    all.t.d2 = c(all.t.d2,t.d2)
    all.mstar = c(all.mstar,rep(length(t.xij),length(t.xij)))
    all.a = c(all.a, cbind(tmp[1:length(t.xij), which(colnames(dat)=="a1")]))

  }

  ugap1 = data.frame(tgtime = all.t.xij, delta = all.t.d1, a1 = all.a, mstar=all.mstar)
  ugap2 = data.frame(tgtime = all.t.zij, delta = all.t.d2, a1 = all.a, mstar=all.mstar)

  #order
  ugap1=ugap1[order(ugap1$tgtime,decreasing=TRUE),]
  ugap2=ugap2[order(ugap2$tgtime,decreasing=TRUE),]
  out=list(n=n,ugap1=ugap1,ugap2=ugap2)

  return(out)
}

##-----point estimation
##rev-biv
RE.biv1=function(beta,dat) {
  mdat=m.dat.chang1(dat,beta)
  n=mdat$n
  ugap1=mdat$ugap1
  ugap2=mdat$ugap2

  ss10=cumsum(1/ugap1$mstar/n)
  ss11=cumsum(ugap1$a1/ugap1$mstar/n)
  sub1=ugap1$delta*(ugap1$a1-ss11/ss10)/ugap1$mstar/sqrt(n)

  ss20=cumsum(1/ugap2$mstar/n)
  ss21=cumsum(ugap2$a1/ugap2$mstar/n)
  sub2=ugap2$delta*(ugap2$a1-ss21/ss20)/ugap2$mstar/sqrt(n)

  out=c(sum(sub1),sum(sub2))
  return(out)
}

RE.uf1=function(beta,dat) {
  tmp.out=RE.biv1(beta,dat)
  out=tmp.out%*%tmp.out
  return(out)
}

RE.uest1=function(init,dat) {
  res=optim(init, RE.uf1, dat=dat, control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}

##-----variance estimation
############################################
#Zeng
############################################
v.est1=function(beta, dat, R)
  #----------------------------------------------------------------------------------------------------------------------------
# first step of variance estimate: estimate V by bootstrap, only need to evaluate the estimating function, no need to solve it
#----------------------------------------------------------------------------------------------------------------------------
{
  id=dat$id
  ids=unique(dat$id)
  n=length(ids)
  freq=table(dat$id)
  index=cumsum( c(0, freq[-n]) )
  p=length(beta)

  A=matrix(rep(NA,R*p),ncol=p)
  for (i in 1:R)
  {
    w=table(sample(ids,n,replace=TRUE))
    s=as.numeric(names(w)) # because of this line, id must be 1:n
    w=as.numeric(w)
    location=NULL
    newid=NULL
    for(ss in 1:length(s))
    {
      location=c(location, rep(index[s[ss]]+(1:freq[s[ss]]),times=w[ss]))
      # since the same subject may be drawn multiple times, new id need to be created to distinguish different duplicates
      # e.g., the first duplicate's id will be original id+1000, the next will be id+2000, etc.
      newid=c(newid,rep(id[index[s[ss]]+(1:freq[s[ss]])],times=w[ss])+rep(((1:w[ss])-1)*1000,each=freq[s[ss]]))
    }
    dat.boot=dat[location,]
    dat.boot$id=newid
    A[i,]=RE.biv1(dat=dat.boot, beta=beta)
  }
  v=cov(A)
  return(v)
}

############################################
#parzen
############################################
#should do v.est first
############################################

RE.bivR1=function(beta,dat,R) {
  mdat=m.dat.chang1(dat,beta)
  n=mdat$n
  # p=length(beta)
  ugap1=mdat$ugap1
  ugap2=mdat$ugap2

  ss10=cumsum(1/ugap1$mstar/n)
  ss11=cumsum(ugap1$a1/ugap1$mstar/n)
  sub1=ugap1$delta*(ugap1$a1-ss11/ss10)/ugap1$mstar/sqrt(n)

  ss20=cumsum(1/ugap2$mstar/n)
  ss21=cumsum(ugap2$a1/ugap2$mstar/n)
  sub2=ugap2$delta*(ugap2$a1-ss21/ss20)/ugap2$mstar/sqrt(n)

  out1=sum(sub1)-R[1]
  out2=sum(sub2)-R[2]
  out=c(out1,out2)
  return(out)
}

RE.ufR1=function(beta,dat,R) {
  tmp.out = RE.bivR1(beta,dat,R)
  out = tmp.out%*%tmp.out
  return(out)
}

RE.uestR1=function(init,dat,R) {
  res=optim(init, RE.ufR1, dat=dat, R=R,control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}

sd.estpar1=function(init,dat,v, B) {
  p=length(init)
  A=matrix(rep(NA,B*p),ncol=p)
  i=0
  while (i < B)
  {
    R=mvrnorm(1,rep(0,p),v)
    est.R=RE.uestR1(init,dat,R)
    if (est.R$conv!=0) next
    i=i+1
    A[i,]=est.R$par
  }
  var_est=cov(A,A) #cov compute the cov between columns
  out=sqrt(diag(var_est))
  return(list(sd=out, covmat=var_est))
}

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for univariate fits using semiparametric regression method on a biv.rec object
#'
#' @description
#' This function fits the model using Chang's Method given one covariate. Called from biv.rec.fit(). No user interface.
#' @param new_data An object that has been reformatted for fit using the biv.rec.reformat() function. Passed from biv.rec.fit().
#' @param cov_names A string with the name of the covariate. Passed from biv.rec.fit().
#' @param CI Passed from biv.rec.fit().
#'
#' @return A list with estimates, SE and variance-covariance matrix.
#'
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats optimize
#' @importFrom stats qnorm
#' @importFrom stats cov
#' @importFrom MASS mvrnorm
#' @importFrom survival Surv
#'
#' @noRd
#' @keywords internal


#one variable regression analysis-Chang's method
chang_univariate <- function(new_data, cov_names, SE) {

  print(paste("Fitting model with covariate", cov_names))
  beta <- rep(0, 2)
  a_index <- which(colnames(new_data)==cov_names)
  colnames(new_data) = c(colnames(new_data)[-a_index], "a1")

  #solve first equation to get all  estimates
  chang1 <- RE.uest1(beta, new_data)

  if (chang1$conv!=0) {
    stop("Max iterations reached. Did not converge.")
  }

  if (is.null(SE)==TRUE) {
    changfit <- data.frame(chang1$par)
    colnames(changfit) <- c("Estimate")
    rownames(changfit) <- c(paste("xij", cov_names), paste("yij", cov_names))
    return(fit=changfit)

  } else {

    print("Estimating standard errors")

    #estimate covariance matrix / std. errors
    chang1v <- v.est1(chang1$par, new_data, R=100)
    chang1sd <- sd.estpar1(beta, new_data, chang1v ,B=50)

    #join all info, put in nice table
    changfit <- data.frame(chang1$par, chang1sd$sd)
    colnames(changfit) <- c("Estimate", "SE")
    rownames(changfit) <- c(paste("xij", cov_names), paste("yij", cov_names))
    return(list(fit=changfit, vcovmat = chang1sd$covmat))
  }

}
