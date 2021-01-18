###########################################################################
############## FUNCTIONS FOR REFERENCE BY MAIN - NOT FOR USER #############
###########################################################################

#             m.dat.chang, all RE, v.est and sd.estpar FUNCTIONS               #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (June, 2018)                          #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________

######################
##-----reformat dataset
m.dat.chang=function(dat,beta) {
  n=length(unique(dat$id))
  mc=max(dat$epi)-1
  p=length(beta)/2
  beta1=beta[1:p]
  beta2=beta[(p+1):(2*p)]
  maxb=apply(cbind(beta1,beta2),1,max)
  #amat=cbind(dat$a1,dat$a2,dat$a3)  changed below
  amat_indexes <- c(9:ncol(dat))
  amat <- as.matrix(dat[,amat_indexes])
  dat$txij=dat$xij*exp(-amat%*%beta1)
  dat$tzij=dat$txij+dat$yij*exp(-amat%*%beta2)
  dat$tci=dat$ci*exp(-amat%*%maxb)

  all.t.xij=all.t.zij=all.t.d1=all.t.d2=all.mstar=all.a=NULL

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
    all.t.xij=c(all.t.xij,t.xij)
    all.t.zij=c(all.t.zij,t.zij)
    all.t.d1=c(all.t.d1,t.d1)
    all.t.d2=c(all.t.d2,t.d2)
    all.mstar=c(all.mstar,rep(length(t.xij),length(t.xij)))
    all.a=rbind(all.a, tmp[1:length(t.xij), amat_indexes])
  }

  ugap1 = cbind(tgtime = all.t.xij, delta = as.integer(all.t.d1), all.a, mstar=all.mstar)
  ugap2 = cbind(tgtime = all.t.zij, delta = as.integer(all.t.d2), all.a, mstar=all.mstar)

  #order
  ugap1=ugap1[order(ugap1$tgtime,decreasing=TRUE),]
  ugap2=ugap2[order(ugap2$tgtime,decreasing=TRUE),]
  out=list(n=n,ugap1=ugap1,ugap2=ugap2)

  return(out)
}

##-----point estimation
##rev-biv
RE.biv=function(beta,dat) {
  mdat=m.dat.chang(dat,beta)
  n=mdat$n
  ugap1=mdat$ugap1
  ugap2=mdat$ugap2
  ncov = length(c(9:ncol(dat)))
  a_indexes1 = a_indexes2 = seq(3, 2+ncov, 1)

  ss10=cumsum(1/ugap1$mstar/n)
  #ss11=apply(cbind(ugap1$a1,ugap1$a2,ugap1$a3)/ugap1$mstar/n,2,cumsum)
  #sub1=ugap1$delta*(cbind(ugap1$a1,ugap1$a2,ugap1$a3)-ss11/ss10)/ugap1$mstar/sqrt(n)
  #changed below
  ss11=apply(ugap1[, a_indexes1]/ugap1$mstar/n,2,cumsum)
  sub1=ugap1$delta*(ugap1[, a_indexes1]-ss11/ss10)/ugap1$mstar/sqrt(n)

  ss20=cumsum(1/ugap2$mstar/n)
  #ss21=apply(cbind(ugap2$a1,ugap2$a2,ugap2$a3)/ugap2$mstar/n,2,cumsum)
  #sub2=ugap2$delta*(cbind(ugap2$a1,ugap2$a2,ugap2$a3)-ss21/ss20)/ugap2$mstar/sqrt(n)
  #changed below
  ss21=apply(ugap2[, a_indexes2]/ugap2$mstar/n,2,cumsum)
  sub2=ugap2$delta*(ugap2[, a_indexes2]-ss21/ss20)/ugap2$mstar/sqrt(n)

  out1=apply(sub1,2,sum)
  out2=apply(sub2,2,sum)
  out=c(out1,out2)
  return(out)
}

RE.uf=function(beta,dat) {
  tmp.out=RE.biv(beta,dat)
  out=tmp.out%*%tmp.out
  return(out)
}

RE.uest=function(init,dat) {
  res=optim(init, RE.uf, dat=dat, control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}


##-----variance estimation
############################################
#Zeng
############################################
v.est=function(beta,dat,R)
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
    #cbind(dat.boot,newid=newid)
    dat.boot$id=newid
    A[i,]=RE.biv(beta,dat.boot)
  }
  v=cov(A)
  return(v)
}

############################################
#parzen
############################################
#should do v.est first
############################################

RE.bivR=function(beta,dat,R) {
  mdat=m.dat.chang(dat,beta)
  n=mdat$n
  p=length(beta)
  ugap1=mdat$ugap1
  ugap2=mdat$ugap2
  ncov = length(c(9:ncol(dat))) #added line
  a_indexes1 = a_indexes2 = seq(3, 2+ncov, 1) #added line

  ss10=cumsum(1/ugap1$mstar/n)
  #ss11=apply(cbind(ugap1$a1,ugap1$a2,ugap1$a3)/ugap1$mstar/n,2,cumsum)
  #sub1=ugap1$delta*(cbind(ugap1$a1,ugap1$a2,ugap1$a3)-ss11/ss10)/ugap1$mstar/sqrt(n)
  ss11=apply(ugap1[, a_indexes1]/ugap1$mstar/n,2,cumsum)
  sub1=ugap1$delta*(ugap1[, a_indexes1]-ss11/ss10)/ugap1$mstar/sqrt(n)

  ss20=cumsum(1/ugap2$mstar/n)
  #ss21=apply(cbind(ugap2$a1,ugap2$a2,ugap2$a3)/ugap2$mstar/n,2,cumsum)
  #sub2=ugap2$delta*(cbind(ugap2$a1,ugap2$a2,ugap2$a3)-ss21/ss20)/ugap2$mstar/sqrt(n)
  #changed below
  ss21=apply(ugap2[, a_indexes2]/ugap2$mstar/n,2,cumsum)
  sub2=ugap2$delta*(ugap2[, a_indexes2]-ss21/ss20)/ugap2$mstar/sqrt(n)

  out1=apply(sub1,2,sum)-R[1:(p/2)]
  out2=apply(sub2,2,sum)-R[(p/2+1):p]
  out=c(out1,out2)
  return(out)
}

RE.ufR=function(beta,dat,R) {
  tmp.out=RE.bivR(beta,dat,R)
  out=tmp.out%*%tmp.out
  return(out)
}

RE.uestR=function(init,dat,R) {
  res=optim(init,RE.ufR,dat=dat,R=R, control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}

sd.estpar=function(init, dat, v, B) {
  p=length(init)
  A=matrix(rep(NA,B*p),ncol=p)
  i=0
  while (i < B)
  {
    R=MASS::mvrnorm(1,rep(0,p),v)
    est.R=RE.uestR(init,dat,R)
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
#' A Function for multivariate fits using semiparametric regression method on a biv.rec object
#'
#' @description
#' This function fits the model using Chang's Method given multiple  covariates. Called from biv.rec.fit(). No user interface.
#' @param new_data An object that has been reformatted for fit using the biv.rec.reformat() function. Passed from biv.rec.fit().
#' @param cov_names A vector with the names of the covariates. Passed from biv.rec.fit().
#' @param SE Passed from biv.rec.fit().
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

#multivariable regression analysis-Chang's method
chang_multivariate <- function(new_data, cov_names, SE) {

  print(paste("Fitting model with covariates:", stringr::str_c(cov_names, collapse = ", "), sep=" "))
  beta <- rep(0, length(cov_names)*2)

  #solve to get all  estimates
  chang <- RE.uest(init = beta, dat=new_data)

  if (chang$conv!=0) {
    stop("Max iterations reached. Did not converge.")
  }

  if (is.null(SE)==TRUE) {
    #return only point estimates
    changfit <- data.frame(chang$par)
    colnames(changfit) <- c("Estimate")
    rownames(changfit) <- c(paste("xij", cov_names), paste("yij", cov_names))
    return(list(fit = as.matrix(changfit)))

  } else {

    print("Estimating standard errors")

    #estimate covariance matrix / std. errors using Parzen's method
    changv <- v.est(chang$par,new_data,R=50)
    changsd <- sd.estpar(init = beta, dat = new_data, v = changv, B=30)

    #Join all info, put in nice table
    changfit <- data.frame(chang$par, changsd$sd)
    colnames(changfit) <- c("Estimate", "SE")
    rownames(changfit) <- c(paste("xij", cov_names), paste("yij", cov_names))
    return(list(fit = as.matrix(changfit), vcovmat = changsd$covmat))

  }

}
