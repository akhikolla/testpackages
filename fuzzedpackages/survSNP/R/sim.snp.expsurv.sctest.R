sim.snp.expsurv.sctest <-
function(n,gtprev,lam,a,b,ztest,diag=FALSE)
  {
    ng=as.integer(rmultinom(1,n,gtprev))
    atime=rexp(n,rep(lam,ng))
    ctime=runif(n,a,b)
    otime=pmin(atime,ctime)
    event=as.integer(atime<ctime)
    SNP=rep(ztest,ng)
    
    mod=summary(coxph(Surv(otime,event)~SNP))
    if(diag)
      print(mod)
    c(event=mean(event),pval=mod$sctest[3])
  }
  
  
