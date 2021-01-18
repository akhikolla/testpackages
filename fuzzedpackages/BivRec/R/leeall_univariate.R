##------symmetric O function
#o.fun=function(t,s,L) {log(min(max(t,s),L))-log(L)}

##-----estimation functions

##proposed method
Pro.ee1=function(beta1, mdat, amat) {
  n=mdat$n
  xmat=mdat$xmat
  delta1=mdat$delta1
  g1mat=mdat$g1mat
  l1=mdat$l1
  mstar=mdat$mstar

  tmp.out=NULL
  for (i in 1:n) {
    A=(amat)-amat[i]
    expA=exp(A*beta1)
    di <- delta1[i,1:mstar[i]]
    xmati <- xmat[i,1:mstar[i]]
    gmati <- g1mat[i,1:mstar[i]]
    subsum <- r2f.pro.ee1(n, nparams=1, di, xmati, gmati, L=l1, expA, subsum = rep(0,n), kcount=mstar[i])
    #subsum=sapply(expA,function(x)mean(delta1[i,1:mstar[i]]*sapply(xmat[i,1:mstar[i]],function(t)o.fun(t,x*t,l1))/g1mat[i,1:mstar[i]]))
    tmp.out=c(tmp.out,sum(A*subsum))
  }
  out=sum(tmp.out)/(n^2)
  return(out)
}

Pro.uf1=function(beta1, mdat, amat) {
  tmp.out=Pro.ee1(beta1,mdat, amat)
  out=tmp.out%*%tmp.out
  return(out)
}

Pro.uest1=function(int, mdat, amat) {
  res=optimize(Pro.uf1,interval=int,mdat=mdat, amat=amat)
  return(list(par=res$minimum,value=res$objective))
}

Pro.ee2=function(beta2, beta1, mdat, amat) {
  n=mdat$n
  xmat=mdat$xmat
  ymat=mdat$ymat
  delta2=mdat$delta2
  g2mat=mdat$g2mat
  l2=mdat$l2
  mstar=mdat$mstar

  tmp.out=NULL
  for (i in 1:n) {
    A=(amat)-amat[i]
    expA1=exp(A*beta1)
    expA2=exp(A*beta2)
    expA=cbind(expA1,expA2)
    di <- delta2[i,1:mstar[i]]
    xmati <- xmat[i,1:mstar[i]]
    ymati <- ymat[i,1:mstar[i]]
    gmati <- g2mat[i,1:mstar[i]]
    subsum <- r2f.pro.ee2(n, nparams=1, di, xmati, ymati, gmati, L=l2, expA, subsum = rep(0,n), kcount=mstar[i])
    # subsum=apply(expA,1,function(x)mean(delta2[i,1:mstar[i]]*apply(cbind(xmat[i,1:mstar[i]],ymat[i,1:mstar[i]]),1,function(t)o.fun(sum(t),x[1]*t[1]+x[2]*t[2],l2))/g2mat[i,1:mstar[i]]))
    tmp.out=c(tmp.out,sum(A*subsum))
  }
  out=sum(tmp.out)/(n^2)
  return(out)
}

Pro.uf2=function(beta2, beta1, mdat, amat) {
  tmp.out=Pro.ee2(beta2, beta1, mdat, amat)
  out=tmp.out%*%tmp.out
  return(out)
}

Pro.uest2=function(int, beta1, mdat, amat) {
  res=optimize(Pro.uf2, interval=int, beta1=beta1, mdat=mdat, amat=amat)
  return(list(par=res$minimum,value=res$objective))
}


##variance estimation
var.est=function(beta1, beta2, mdat, amat) {
  n=mdat$n
  mc=mdat$mc
  xmat=mdat$xmat
  ymat=mdat$ymat
  delta1=mdat$delta1
  delta2=mdat$delta2
  g1mat=mdat$g1mat
  g2mat=mdat$g2mat
  l1=mdat$l1
  l2=mdat$l2
  mstar=mdat$mstar

  xi=matrix(0,length(c(beta1,beta2)),length(c(beta1,beta2)))
  gam1=gam21=gam22=rep(0,length(beta1))
  for (i in 1:n) {
    A=(amat)-amat[i]
    expA1=exp(A*beta1)
    expA2=exp(A*beta2)
    expA=cbind(expA1,expA2)

    d1i <- delta1[i,1:mstar[i]]
    d2i <- delta2[i,1:mstar[i]]
    xmati <- xmat[i,1:mstar[i]]
    ymati <- ymat[i,1:mstar[i]]
    gmati1 <- g1mat[i,1:mstar[i]]
    gmati2 <- g2mat[i,1:mstar[i]]

    subsum <- rep(0,n)
    sub1.xi1 <- r2f.pro.ee1(n, nparams=1, di=d1i, xmati, gmati=gmati1, L=l1, expA=expA1, subsum, kcount=mstar[i])
    sub1.xi2 <- r2f.pro.ee2(n, nparams=1, di=d2i, xmati, ymati, gmati=gmati2, L=l2, expA, subsum, kcount=mstar[i])

    #sub1.xi1=sapply(expA1,function(x) mean(delta1[i,1:mstar[i]]*sapply(xmat[i,1:mstar[i]],function(t)o.fun(t,x*t,l1))/g1mat[i,1:mstar[i]]))
    #sub1.xi2=apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*apply(cbind(xmat[i,1:mstar[i]],ymat[i,1:mstar[i]]),1,function(t)o.fun(sum(t),x[1]*t[1]+x[2]*t[2],l2))/g2mat[i,1:mstar[i]]))

    sub2 <- r2f.pro.var(n, nparams=1, xmat, ymat, gmatx=g1mat, gmaty=g2mat, l1, l2,
                         expAx=expA1, expAy=expA2, subsumx=subsum, subsumy=subsum, dx=delta1, dy=delta2, mstar, mc)
    sub2.xi1 <- sub2[,1]
    sub2.xi2 <- sub2[,2]

    #sub2.xi1=sub2.xi2=rep(0,n)
    # for (j in 1:n) {
    #   sub2.xi1[j]=mean(delta1[j,1:mstar[j]]*sapply(xmat[j,1:mstar[j]],function(t)o.fun(t,t/expA1[j],l1))/g1mat[j,1:mstar[j]])
    # }
    # for (j in 1:n) {
    #   sub2.xi2[j]=mean(delta2[j,1:mstar[j]]*apply(cbind(xmat[j,1:mstar[j]],ymat[j,1:mstar[j]]),1,function(t)o.fun(sum(t),t[1]/expA1[j]+t[2]/expA2[j],l2))/g2mat[j,1:mstar[j]])
    # }

    tmp.xi1=sum(A*(sub1.xi1-sub2.xi1))/(n^(3/2))
    tmp.xi2=sum(A*(sub1.xi2-sub2.xi2))/(n^(3/2))

    xi=xi+c(tmp.xi1,tmp.xi2)%o%c(tmp.xi1,tmp.xi2)

    Amat=sapply(A,function(x) x%o%x)
    tmp.sub.gam1=apply(cbind(xmat[i,1],expA1*xmat[i,1],l1),1,function(x)(x[1]<=x[2])*(max(x[1],x[2])<=x[3]))
    tmp.sub.gam2=apply(cbind(xmat[i,1]+ymat[i,1],expA1*xmat[i,1]+expA2*ymat[i,1],l2),1,function(x)(x[1]<=x[2])*(max(x[1],x[2])<=x[3]))
    sub.gam1=(Amat)*tmp.sub.gam1*mean(delta1[i,1:mstar[i]]/g1mat[i,1:mstar[i]])
    sub.gam21=(Amat)*tmp.sub.gam2*apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*(x[1]*xmat[i,1:mstar[i]])/((x[1]*xmat[i,1:mstar[i]]+x[2]*ymat[i,1:mstar[i]])*g2mat[i,1:mstar[i]])))
    sub.gam22=(Amat)*tmp.sub.gam2*apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*(x[2]*ymat[i,1:mstar[i]])/((x[1]*xmat[i,1:mstar[i]]+x[2]*ymat[i,1:mstar[i]])*g2mat[i,1:mstar[i]])))

    gam1=gam1+sum(sub.gam1)/(n^2)
    gam21=gam21+sum(sub.gam21)/(n^2)
    gam22=gam22+sum(sub.gam22)/(n^2)
  }
  gam1=matrix(gam1,length(beta1),length(beta1))
  gam21=matrix(gam21,length(beta2),length(beta2))
  gam22=matrix(gam22,length(beta2),length(beta2))
  gamm=rbind(cbind(gam1,matrix(0,length(beta1),length(beta2))),cbind(gam21,gam22))

  mat=solve(gamm)%*%xi%*%t(solve(gamm))

  se1=sqrt(diag(mat)/n)[1]
  se2=sqrt(diag(mat)/n)[2]
  return(list(se1=se1,se2=se2, vcovmat=mat/n))
}

##################### FUNCTION NOT FOR USER #######################
###################################################################
#' A Function for univariate fits using semiparametric regression method on a bivrecSurv object
#'
#' @description
#' This function fits the semiparametric model given one  covariate. Called from bivrecReg(). No user interface.
#' @param response Passed from bivrecReg().
#' @param amat Passed from bivrecReg().
#' @param cov_names Passed from bivrecReg().
#' @param SE Passed from bivrecReg()
#' @return A list with estimates, SE and variance-covariance matrix.
#'
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats optimize
#' @importFrom stats qnorm
#' @importFrom stringr str_c
#'
#' @noRd
#' @keywords internal

#MAIN PROGRAM FOR univariate regression analysis

#multivariable regression analysis
leeall_univariate <- function(response, amat, cov_names, SE){

  print(paste("Fitting model with covariate", cov_names))

  #solve first equation to get beta1 values - related to xij
  pro1 <- Pro.uest1(c(-2,2), mdat=response, amat=amat)[[1]]

  #solve second equation to get beta2 values - related to yij
  pro2 <- Pro.uest2(c(-2,2), pro1, mdat=response, amat=amat)[[1]]

  if (SE==TRUE) {

    print("Estimating standard errors")
    #estimate covariance matrix and get diagonal then std. errors
    sd_est=var.est(pro1, pro2, mdat=response, amat = amat)
    univ_fits <- data.frame(c(pro1, pro2), c(sd_est[[1]],sd_est[[2]]))
    colnames(univ_fits) <- c("Estimate", "SE")
    rownames(univ_fits) <- c(paste("xij", cov_names), paste("yij", cov_names))
    result <- list(fit = as.matrix(univ_fits),  vcovmat = sd_est[[3]])

  } else {
    #return point estimates only
    univ_fits <- data.frame(c(pro1, pro2))
    colnames(univ_fits) <- c("Estimate")
    rownames(univ_fits) <- c(paste("xij", cov_names), paste("yij", cov_names))
    result <- list(pro1$par, pro2$par, univ_fits)

  }

  return(result)
}
