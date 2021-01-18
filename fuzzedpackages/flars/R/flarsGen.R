GQ=0

FccaXYdir <- function(Lx, Lp, LPhi, resp, l1, l2, cv=1) {
  .Call('flars_FccaXYdir',Lx, Lp, LPhi, resp, l1, l2, cv=1, PACKAGE = 'flars')
}

FccaXYdir0 <- function(Lx, Lp, LPhi, resp, l1, l2, cv=1) {
  .Call('flars_FccaXYdir0',Lx, Lp, LPhi, resp, l1, l2, cv=1, PACKAGE = 'flars')
}
cov.linear=function(hyper,Data,Data.new=NULL){
  if(is.null(Data.new)) Data.new=Data
  hyper=lapply(hyper,exp);n.hyper=length(hyper$linear.a)
  cov.lin=xixj(Data.new,Data,hyper$linear.a)
  return(cov.lin)
}

cov.pow.ex=function(hyper,Data,Data.new=NULL,gamma=1){
  #hyper is a list of hyper-parameters
  if(is.null(gamma)) gamma=1
  hyper=lapply(hyper,exp);
  datadim=dim(Data)

  v.power=xixj_sta(Data,Data.new,hyper$pow.ex.w,power=gamma)
  
  exp.v.power=hyper$pow.ex.v*exp(-v.power/2)
}

xixj=function(mat,mat.new=NULL,a=NULL){
  mat=as.matrix(mat)
  mdim=dim(mat)

  if(is.null(mat.new)){
    mat.new=mat
  }
  
  if(is.null(a))  a=rep(1,mdim[2])
  if(length(a)<mdim[2]) {
    a1=rep(1,mdim[2])
    a1[1:length(a)]=a
    a=a1;rm(a1)
    warning('number of "a" is less than the number of columns, use 1 as the missing "a"')
  }
  if(length(a)>mdim[2]) {
    a=a[1:mdim[2]]
    warning('number of "a" is more than the number of columns, omit the extra "a"')
  }
  
  aa=matrix(rep(a,mdim[1]),ncol=mdim[2],byrow=T)
  out=(aa*mat)%*%t(mat.new)
  return(out)
}

xixj_sta=function(mat,mat.new=NULL,w=NULL,power=NULL){
  mat=as.matrix(mat)
  if(is.null(mat.new)) mat.new=mat
  mdim=dim(mat);mdim.new=dim(mat.new)
  cov.=matrix(sapply(1:mdim[1],function(i) matrix(rep(mat[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-mat.new),ncol=mdim[1])
  
  if(is.null(power)) power=1
  cov.=((cov.)^2)^power;
  if(is.null(w)) {
    w=rep(1,mdim[2])
    warning('missing "weight", use 1 instead')
  }
  if(length(w)==1&mdim[2]>1){
    w=rep(w,mdim[2])
    warning('only one "weight" found, applied to all columns')
  }
  
  if(length(w)>1&length(w)<mdim[2]){
    w1=rep(1,mdim[2])
    w1[1:length(w)]=w
    w=w1;rm(w1)
    warning('number of "weight" is less than the number of columns, use 1 as the missing "weight"')
  }
  if(length(w)>mdim[2]){
    w=w[1:mdim[2]]
    warning('number of "weight" is more than the number of columns, omit the extra "weight"')
  }
  
  wmat=matrix(rep(w,each=dim(cov.)[1]*dim(cov.)[2]/mdim[2]),ncol=dim(cov.)[2],byrow=T)
  
  cov.=wmat*cov.
  
  cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
  if(mdim[2]>1){
    for(i in 1:(mdim[2]-1)){
      cov..=cov..+cov.[1:mdim.new[1],,drop=F];cov.=cov.[-(1:mdim.new[1]),,drop=F]}
    cov.=cov..+cov.
  }
  return(cov.)  
}



ind_fv_gen=function(effect_sample_size=10,num_fv=3,nsample=80,seed=NULL){
  ## generate independent functional variables
  if(is.null(seed)) seed=sample(seq(1,1e4),size=1)
  set.seed(seed)
  
  x=matrix(rnorm(nsample*num_fv*effect_sample_size),ncol=effect_sample_size*num_fv)
  idx_factor=rep(rep(1:num_fv,each=nsample),effect_sample_size)
  x=lapply(split(x,idx_factor),function(i) matrix(i,ncol=effect_sample_size))
  
  for(i in 2:num_fv){
    a=do.call('cbind',x[1:(i-1)])
    b=x[[i]]
    out=resid(lm(b~a))
    x[[i]]=out
  }
  
  return(x)
}

mat2fd=function(mat,fdList=NULL){
  fl=list(time=seq(0,1,len=ncol(mat)),nbasis=min(as.integer(ncol(mat)/5),23),norder=6,bSpline=TRUE,Pen=c(0,0),lambda=1e-4)
  nbasis=c(fdList$nbasis,fl$nbasis)[1]
  norder=c(fdList$norder,fl$norder)[1]
  lambda=c(fdList$lambda,fl$norder)[1]
  bSpline=c(fdList$bSpline,fl$bSpline)[1]
  time=list(a=fdList$time,b=fl$time)
  time=time[[which(unlist(lapply(time,is.null))^2==0)[1]]]
  if(1-bSpline) fl$Pen=c(c(0,(2*pi/diff(range(time)))^2,0))
  Pen=list(a=fdList$Pen,b=fl$Pen)
  Pen=Pen[[which(unlist(lapply(Pen,is.null))^2==0)[1]]]
  if(bSpline)  basis=create.bspline.basis(range(time),nbasis,norder)
  if(1-bSpline) basis=create.fourier.basis(range(time),nbasis,diff(range(time)))
  Par=vec2Lfd(Pen,range(time))
  matfd=smooth.basisPar(time,t(mat),basis,Lfdobj=Par,lambda)$fd
  return(matfd)
}

add_shape=function(mat,vec){
  dim_vec=length(vec)
  if(nrow(mat)==dim_vec & nrow(mat)!=ncol(mat)) mat=t(mat)
  mat=apply(mat,1,function(iii) return(iii+vec))
  if(nrow(mat)==dim_vec & nrow(mat)!=ncol(mat)) mat=t(mat)
  return(mat)
}



f_var_from_f_independent=function(f_independent,weight,intercept=NULL, slope=NULL){
  fi=f_independent
  nf=length(fi)
  ni=ns=0
  if(!is.null(intercept)) ni=length(intercept)
  if(!is.null(slope)) ns=length(slope)
  
  if(is.null(intercept)) intercept=0
  if(is.null(slope)) slope=0
  
  z=vector('list',length=nrow(weight))
  weight=split(weight,rep(1:nrow(weight),ncol(weight)))
  for(i in seq_along(z)){
    z0=lapply(1:nf,function(k)
      fi[[k]]*weight[[i]][k]*matrix(slope[[sample(ns,1)]],ncol=ncol(fi[[k]]),nrow=nrow(fi[[k]]),byrow = T)+matrix(intercept[[sample(ni,1)]],ncol=ncol(fi[[k]]),nrow=nrow(fi[[k]]),byrow = T))
    z[[i]]=Reduce('+',z0)
  }
  return(z)
}

rmse=function(t,a){ 
  #compute the root mean squar error between two vectors
  y = sqrt(mean((a-t)^2,na.rm = T))
  return(y)
}

getL=function(dimension,order){
  l=matrix(0,nrow=dimension-order,ncol=dimension)
  for(i in 1:nrow(l))
    l[i,i+0:order]=(c(-1,1,0)+(order-1)*c(0,1,-1))[1:(order+1)]
  return(l)
}

get2L=function(dimension,DelT){
  l=matrix(0,nrow = dimension-2,ncol=dimension)
  for(i in 1:nrow(l)){
    l[i,i+0:2]=2/c((DelT[i]+DelT[i+1])*DelT[i],-DelT[i]*DelT[i+1],(DelT[i]+DelT[i+1])*DelT[i+1])
  }
  return(l)
}

cp=function(mat1,mat2=NULL){
  if(is.data.frame(mat1)) mat1=data.matrix(mat1)
  if(is.data.frame(mat2)) mat2=data.matrix(mat2)
  return(crossprod(mat1,mat2))
}
tcp=function(mat1,mat2=NULL){
  if(is.data.frame(mat1)) mat1=data.matrix(mat1)
  if(is.data.frame(mat2)) mat2=data.matrix(mat2)
  return(tcrossprod(mat1,mat2))
}

GenWeight=function(nx,seed,avoid){
  a=rnorm(nx)
  a[avoid]=rnorm(length(avoid),sd = 1e-1)
  a=nx*a/sum(abs(a))
  return(abs(a))
}

plot_fm=function(f_var,dim=NULL){
  nfm=length(f_var)
  if(is.null(dim))  par(mfrow=c(((nfm/2)%%1!=0)+as.integer(nfm/2),2))
  if(!is.null(dim)) par(mfrow=c(dim[1],dim[2]))
  for(i in 1:nfm){
    fv=f_var[[i]]
    plot(-100,-100,xlim=c(0,ncol(fv)),ylim=range(fv),
         main=substitute( paste("x"[{x000}],"(t)",sep=''), list(x000=i) ),
         ylab='',xlab='index')
    for(ii in 1:nrow(fv)) lines(1:ncol(fv),fv[ii,],col=col1(1.2*nrow(fv))[ii])
  }
  par(mfrow=c(1,1))
}

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", 
                           "#007FFF", "blue", "#00007F"))




ran_fv_gen=function(effect_sample_size=10,num_fv=3,nsample=80,seed=NULL){
  ## generate independent functional variables
  if(is.null(seed)) seed=sample(seq(1,1e4),size=1)
  set.seed(seed)
  
  x=matrix(rnorm(nsample*num_fv*effect_sample_size),ncol=effect_sample_size*num_fv)
  idx_factor=rep(rep(1:num_fv,each=nsample),effect_sample_size)
  x=lapply(split(x,idx_factor),function(i) matrix(i,ncol=effect_sample_size))
  
  
  
  return(x)
}

data_generation=function(seed,nsamples=80,hyper=NULL,var_type=c('f','m'),cor_type=1:6,uncorr=TRUE,nVar=8){
  set.seed(seed)
  
  var_type=var_type[1]
  cor_type=cor_type[1]
  
  data.gen=T
  
  if(is.null(hyper)) 
    hp <- list('pow.ex.w'=log(50),'linear.a'=log(2),'pow.ex.v'=log(20),'vv'=log(0.5))
  if(!is.null(hyper)) hp=hyper
  
  if(data.gen){
    ## generate independent random vectors
    # seed=30
    if(uncorr==T)
      x=ind_fv_gen(effect_sample_size=10,num_fv=nVar,nsample=nsamples,seed=seed)
    if(uncorr==F)
      x=ran_fv_gen(effect_sample_size=10,num_fv=nVar,nsample=nsamples,seed=seed)
    effect_sample_size=10;num_fv=nVar;nsample=nsamples
    G=cov.pow.ex(hp,seq(0,1,len=10))+cov.linear(hp,seq(0,1,len=10))
    G2=chol(G)
    x=lapply(x,'%*%',G2)
    
    ## trun to functional case
    f=lapply(x,function(i){
      mat2fd(i,fdList=list(nbasis=2*effect_sample_size,norder=min(effect_sample_size-1,6),lambda=1e-7))
    })
    
    ## discrete into 100 time points
    ntime=100
    fm0=lapply(f,function(i) t(eval.fd(seq(0,1,len=ntime),i)))
    
    # add noise
    err_range=diff(apply(sapply(fm0,range),1,mean))/20
    fm0=lapply(fm0,function(i) i+eval.fd(seq(0,1,len=ntime),mat2fd(t(rnorm(50,sd=err_range)),list(lambda=1e-3)))[,1])
    
    ## define shape
    shape0=list(s1=sin(seq(0,pi,len=ntime)),
                s2=pnorm(seq(0.6,1.2,len=ntime),mean=.8,sd=0.01)+0.2*exp(seq(-1,2.5,len=ntime)),
                s3=cos(seq(0,1.5,len=ntime))^2,
                s4=sin(cos(seq(0,3,len=ntime))),
                s5=cos(sin(seq(-2,3,len=ntime))),
                s6=seq(-2,3,len=ntime)^3-cos(seq(-2,3,len=ntime)^2)^3,
                s7=cos(sqrt(exp(seq(-2,3,len=ntime)))),
                s8=sin(sqrt(exp(seq(-2,3,len=ntime)))))
    
    shape=lapply(shape0,function(i)return(i*err_range*nsample/diff(range(i))))
    
    # add shape
    
    fm=fm0
    for(i in seq_along(fm)) fm[[i]]=add_shape(fm[[i]],shape[[sample(seq_along(shape0),1)]])
    
    b1=princomp(x[[1]], cor=F)$loadings[,1]
    b1=mat2fd(t(as.matrix(b1)),fdList=list(nbasis=5,norder=4,lambda=0))
    bm=eval.fd(seq(0,1,len=ntime),b1)
    b1const=(max(fm[[1]]%*%bm)-min(fm[[1]]%*%bm))
    y01=fm[[1]]%*%bm/b1const
    y1=y01
    
    
    b2=princomp(x[[2]], cor=F)$loadings[,1]
    p=dnorm((rnorm(2*ntime)))
    idx=sort(sample(1:(2*ntime),10,prob=p/sum(p)))
    time2=seq(0,1,len=2*ntime)[idx]
    time2=(time2-(min(time2)))/max(time2)
    b2=mat2fd(t(as.matrix(b2)),fdList=list(time=time2,nbasis=7,norder=4,lambda=0))
    bm=eval.fd(seq(min(time2),max(time2),len=ntime),b2)
    b2const=(max(fm[[2]]%*%bm)-min(fm[[2]]%*%bm))
    y02=fm[[2]]%*%bm/b2const
    y2=y02
    
    noise=rnorm(nsample,mean=0,sd=0.5)
    b3=princomp(x[[3]], cor=F)$loadings[,1]
    b3=mat2fd(t(as.matrix(b3)),fdList=list(nbasis=5,norder=4,lambda=0))
    bm=eval.fd(seq(0,1,len=ntime),b3)
    b3const=(max(fm[[3]]%*%bm)-min(fm[[3]]%*%bm))
    y03=fm[[3]]%*%bm/b3const
    y3=y03
    
    mu=sample(c(0,5,10,20,100),size = 1)
    noise=mu+noise*0.1
  }
  BetaT=list(b1=b1,b2=b2,b3=b3)
  bConst=c(b1const=b1const,b2const=b2const,b3const=b3const)
  
  if(var_type=='f'){
    if(cor_type==1){
      fsm=fm
      y=y1+y2+y3+noise
    }
    if(cor_type==2){
      fmm=list('vector',8)
      fmm[c(1:3,5:8)]=fm[c(1:3,5:8)]
      
      GenWeight=function(nx,seed,avoid){
        a=rnorm(nx)
        a[avoid]=rnorm(length(avoid),sd = 1e-1)
        a=nx*a/sum(abs(a))
        return(abs(a))
      }
      w1=w2=w3=rep(0,8)
      w1[c(1,4)]=GenWeight(2,2,c(1,2))
      w=t(as.matrix(w1))
      z=f_var_from_f_independent(f_independent=fm,weight=w,intercept = shape0,slope = shape0)
      fmm[c(4)]=z
      fsm=fmm
      
      y=y1+y2+y3+noise
      
    }
    
    
    if(cor_type==3){
      fmm=list('vector',8)
      fmm[c(1:3,7,8)]=fm[c(1:3,7,8)]    
      
      GenWeight=function(nx,seed,avoid){
        a=rnorm(nx)
        a[avoid]=rnorm(length(avoid),sd = 1e-1)
        a=nx*a/sum(abs(a))
        return(abs(a))
      }
      w1=w2=w3=rep(0,8)
      w1[c(1,4)]=GenWeight(2,2,c(1,2))
      w2[c(2,5)]=GenWeight(2,3,c(1,2))
      w3[c(3,6)]=GenWeight(2,4,c(1,2))
      w=rbind(w1,w2,w3)
      z=f_var_from_f_independent(f_independent=fm,weight=w,intercept = shape0,slope = shape0)
      fmm[c(4:6)]=z
      fsm=fmm
      
      y=y1+y2+y3+noise
      
    }
    
    
    if(cor_type==4){
      fmm=list('vector',8)
      fmm[c(1:3,6:8)]=fm[c(1:3,6:8)]
      
      GenWeight=function(nx,seed,avoid){
        a=rnorm(nx)
        a[avoid]=rnorm(length(avoid),sd = 1e-1)
        a=nx*a/sum(abs(a))
        return(abs(a))
      }
      w1=w2=rep(0,8)
      w1[c(1,4)]=c(0.2,0.8)
      w2[c(1,5)]=c(0.2,0.8)
      w=rbind(w1,w2)
      z=f_var_from_f_independent(f_independent=fm,weight=w,intercept = shape0,slope = shape0)
      fmm[c(4,5)]=z
      fsm=fmm
      
      y=y1+y2+y3+noise
      
    }
    beta_xs=NULL
    
  }
  
  if(var_type=='m'){
    xs=x[-8]
    xsMat=matrix(0,ncol=5,nrow=nsamples)
    xsMat=apply(xsMat,2,rnorm,n=nsamples)
    if(uncorr==F){
      a=do.call(cbind,xs)
      b=xsMat[,1]
      xsMat[,1]=resid(lm(b~a))
      for(i in 2:ncol(xsMat)){
        a=cbind(do.call('cbind',xs),xsMat[,1:(i-1)])
        b=xsMat[,i]
        out=resid(lm(b~a))
        xsMat[,i]=out
      }
    }
    xsMat=scale(xsMat)
    
    beta_xs=1/apply(apply(xsMat,2,range),2,diff)
    beta_xs[4:5]=0
    y_xs=xsMat%*%beta_xs
    xsList=lapply(split(xsMat,rep(1:ncol(xsMat),each=nrow(xsMat))),as.matrix)
    y=y1+y2+y3+y_xs
    
    if(cor_type==1){
      fsm=fm
      fsm=fsm[-8]
      
      fsm=c(fsm,xsList)
      y=y+noise
    }
    
    if(cor_type==2){
      fmm=list('vector',8)
      fmm[c(1:8)]=fm[1:8]  
      
      GenWeight=function(nx,seed,avoid){
        a=rnorm(nx)
        a[avoid]=rnorm(length(avoid),sd = 1e-1)
        a=nx*a/sum(abs(a))
        return(abs(a))
      }
      w1=w2=rep(0,8)
      w1[c(1,7,5)]=c(0.2,-0.01,0.3)
      w2[c(1,7,4)]=c(0.2,-0.01,0.3)
      w=rbind(w1,w2)
      z=f_var_from_f_independent(f_independent=fm,weight=w,intercept = shape0,slope = shape0)
      fmm[c(4,5)]=z
      
      fsm=fmm
      fsm=fsm[-8]
      
      fsm=c(fsm,xsList)
      y=y+noise
    }
    
    if(cor_type==3){
      fsm=fm
      fsm=fsm[-8]
      
      fsm=c(fsm,xsList)
      fsm[[11]]=0.2*xsList[[1]]+0.8*xsList[[4]]
      fsm[[12]]=0.2*xsList[[1]]+0.8*xsList[[5]]
      
      y=y+noise
    }
    
    if(cor_type==4){
      fmm=list('vector',8)
      fmm[c(1:8)]=fm[1:8]  
      
      GenWeight=function(nx,seed,avoid){
        a=rnorm(nx)
        a[avoid]=rnorm(length(avoid),sd = 1e-1)
        a=nx*a/sum(abs(a))
        return(abs(a))
      }
      w1=w2=rep(0,8)
      w1[c(1,7,5)]=c(0.2,-0.01,0.3)
      w2[c(1,7,4)]=c(0.2,-0.01,0.3)
      w=rbind(w1,w2)
      z=f_var_from_f_independent(f_independent=fm,weight=w,intercept = shape0,slope = shape0)
      fmm[c(4,5)]=z
      
      fsm=fmm
      fsm=fsm[-8]
      
      fsm=c(fsm,xsList)
      fsm[[11]]=0.2*xsList[[1]]+0.8*xsList[[4]]
      fsm[[12]]=0.2*xsList[[1]]+0.8*xsList[[5]]
      
      y=y+noise
    }
    
    if(cor_type==5){
      fmm=list('vector',8)
      fmm[c(1:8)]=fm[1:8]  
      
      GenWeight=function(nx,seed,avoid){
        a=rnorm(nx)
        a[avoid]=rnorm(length(avoid),sd = 1e-1)
        a=nx*a/sum(abs(a))
        return(abs(a))
      }
      w1=w2=rep(0,8)
      w1[c(1,7,5)]=c(0.2,-0.01,0.3)
      w2[c(1,7,4)]=c(0.2,-0.01,0.3)
      w=rbind(w1,w2)
      z=f_var_from_f_independent(f_independent=fm,weight=w,intercept = shape0,slope = shape0)
      fmm[c(4,5)]=z
      
      fsm=fmm
      fsm=fsm[-8]
      
      fsm=c(fsm,xsList)
      fsm[[11]]=0.2*xsList[[1]]+0.8*xsList[[4]]
      fsm[[12]]=0.2*xsList[[1]]+0.8*xsList[[5]]
      
      fsm[[6]]=apply(fsm[[6]],2,function(i6) i6+xsList[[1]])
      
      y=y+noise
    }
    
    if(cor_type==6){
      fmm=list('vector',8)
      fmm[c(1:8)]=fm[1:8]  
      
      GenWeight=function(nx,seed,avoid){
        a=rnorm(nx)
        a[avoid]=rnorm(length(avoid),sd = 1e-1)
        a=nx*a/sum(abs(a))
        return(abs(a))
      }
      w1=w2=rep(0,8)
      w1[c(1,2,8,4)]=c(0.8,0,1.2,0.5)
      w=rbind(w1,w2)
      z=f_var_from_f_independent(f_independent=fm,weight=w,intercept = shape0,slope = shape0)
      fmm[c(4,5)]=z
      
      fsm=fmm
      fsm=fsm[-8]
      
      fsm=c(fsm,xsList)
      fsm[[11]]=0.2*xsList[[1]]+0.8*xsList[[4]]
      fsm[[12]]=0.2*xsList[[1]]+0.8*xsList[[5]]
      
      fsm[[6]]=apply(fsm[[6]],2,function(i6) i6+0.2*xsList[[1]])
      
      y=y+noise
    }  
  }  
  BetaT$S=beta_xs
  
  return(list(x=fsm,y=y,BetaT=BetaT,bConst=bConst,noise=noise,mu=mu))
  
}



fccaGen=function(xL,yVec,type=c('dir','cor','a','all'),method=c('basis','gq','raw'),GCV=TRUE,control=list()){
  if(!is.list(xL)) xL=list(xL)
  
  if(is.list(xL)) xL=lapply(xL,as.matrix)
  
  nList=length(xL)
  ndimL=sapply(xL,ncol)
  idxFD=seq_along(xL)[ndimL>1]
  
  if(length(idxFD)==0){
    x=sapply(xL,as.matrix)
    type=type[1]
    xc=crossprod(x)
    xy=t(x)%*%yVec
    Beta=as.matrix(rep(NA,length(yVec)))
    try(Beta<-solve(xc,xy),T)
    
    if(type=='dir' | type=='a'){
      return(Beta)
    }
    if(type=='cor')
      return(cor(x%*%Beta,yVec))
    if(type=='all'){
      return(list('corr'=cor(x%*%Beta,yVec),'a'=Beta,'PureScalar'=1,'TraceHat'=dim(x)))
    }
  }
  
  
  method=method[1]
  gqWeight=NULL
  
  
  pL=phiL=vector('list',length=nList)
  
  if(method=='basis'){
    con=list(nbasis=18,norder=6,pen1=10^(seq((-20),5,len=41)),pen2=0.01,t=seq(0,1,len=max(sapply(xL,ncol),na.rm = T)))
    con[(namc <- names(control))] <- control
    
    nbasis=con$nbasis;norder=con$norder;pen1=con$pen1
    pen2=con$pen2;t=con$t
    
    if(!is.list(t)){
      tRange=range(t)
      t=lapply(xL,function(i){
        if(ncol(i)==length(t)){
          return(t)
        }
        if(ncol(i)!=length(t) & ncol(i)>1){
          return(seq(tRange[1],tRange[2],len=ncol(i)))
        }
        if(ncol(i)==1){
          return((tRange[1]+tRange[2])/2)
        }
      })
    }
    if(is.list(t)){
      if(length(t)!=nList)
        t=lapply(xL,function(i)return(t[[1]]))
    }
    
    
    if(length(nbasis)!=length(norder)){
      nbasis=nbasis[1]
      norder=norder[1]
      warning('length of nbasis is different from length of norder')
    }
    if(length(nbasis)!=nList)
      nbasis=rep(nbasis[1],nList)
    if(length(norder)!=nList)
      norder=rep(norder[1],nList)
    
    for(iL in idxFD){
      tRange=range(t[[iL]])
      nb=nbasis[iL];no=norder[iL]
      tiL=t[[iL]]
      spline=create.bspline.basis(tRange,nbasis = nb,norder = no)
      MS=eval.basis(tiL,spline,0)
      MS2=eval.basis(tiL,spline,2)
      
      pL[[iL]]=crossprod(MS2)
      phiL[[iL]]=MS/(ncol(xL[[iL]]))
    }
  }
  
  if(method=='gq'){
    con=list(nP=18,pen1=10^(seq((-20),5,len=21)),pen2=0.01,t=seq(-1,1,len=max(sapply(xL,ncol),na.rm = T)))
    con[(namc <- names(control))] <- control
    
    t=con$t;nP=con$nP;pen1=con$pen1;pen2=con$pen2
    
    if(!is.list(t)){
      tRange=range(t)
      t=lapply(xL,function(i){
        if(ncol(i)==length(t)){
          return(t)
        }
        if(ncol(i)!=length(t) & ncol(i)>1){
          return(seq(tRange[1],tRange[2],len=ncol(i)))
        }
        if(ncol(i)==1){
          return((tRange[1]+tRange[2])/2)
        }
      })
    }
    if(is.list(t)){
      if(length(t)!=nList)
        t=lapply(xL,function(i)return(t[[1]]))
    }
    if(GQ==0){
      .GlobalEnv$modelChronic  <-  readRDS(
        file.path(system.file('inst/extdata',package='flars'), 'GQ.rds'))
    }
    
    gqWeight=data.matrix(GQ[GQ$n==nP,c(1,2)])
    gqWeight2=data.matrix(GQ[GQ$n==(nP-2),c(1,2)])
    xi=as.numeric(gqWeight[,1])
    wi=as.matrix(gqWeight[,2])
    wi=wi[order(xi)]
    xi=sort(xi)
    xi2=as.numeric(gqWeight2[,1])
    wi2=as.matrix(gqWeight2[,2])
    wi2=wi2[order(xi2)]
    xi2=sort(xi2)
    idxL=vector('list',length=length(idxFD))
    
    xL0=xL
    for(iL in idxFD){
      dimFx=dim(xL[[iL]])
      ncols=dimFx[2]
      nrows=dimFx[1]
      zi=t[[iL]]
      idx=as.numeric(sapply(xi,function(i) which.min(abs(zi-i))))
      idxL[[iL]]=idx
      xGQ=xL[[iL]][,idx]
      xL[[iL]]=xGQ
      
      pL[[iL]]=t(get2L(ncol(xL[[iL]]),diff(xi)))%*%diag(wi2)%*%get2L(ncol(xL[[iL]]),diff(xi))
      phiL[[iL]]=diag(wi)
    }
  }
  
  if(method=='raw'){
    con=list(pen1=10^(seq((-10),-5,len=21)),pen2=0.01,t=seq(0,1,len=max(sapply(xL,ncol),na.rm = T)))
    con[(namc <- names(control))] <- control
    
    t=con$t;pen1=con$pen1;pen2=con$pen2
    
    if(!is.list(t)){
      tRange=range(t)
      t=lapply(xL,function(i){
        if(ncol(i)==length(t)){
          return(t)
        }
        if(ncol(i)!=length(t) & ncol(i)>1){
          return(seq(tRange[1],tRange[2],len=ncol(i)))
        }
        if(ncol(i)==1){
          return((tRange[1]+tRange[2])/2)
        }
      })
    }
    if(is.list(t)){
      if(length(t)!=nList)
        t=lapply(xL,function(i)return(t[[1]]))
    }
    
    for(iL in idxFD){
      pL[[iL]]=crossprod(get2L(ncol(xL[[iL]]),diff(t[[iL]])))
      phiL[[iL]]=diag(1,ncol(xL[[iL]]))/(ncol(xL[[iL]]))
    }
  }
  
  if(length(idxFD)<nList){
    pL[-idxFD]=lapply(pL[-idxFD],function(i) i=matrix(0))
    phiL[-idxFD]=lapply(phiL[-idxFD],function(i) i=matrix(1))
  }
  
  type=type[1]
  GCV=1
  if(type!='all')
    B=FccaXYdir(xL,pL,phiL,yVec,pen1,pen2,cv=GCV)
  if(type=='all'){
    B_all=FccaXYdir0(xL,pL,phiL,yVec,pen1,pen2,cv=GCV)
    B=B_all$betahat
    B_K=B_all$K
    B_S=B_all$S
    B_lam1=B_all$lambda1
    B_lam2=B_all$lambda2
    GCV_mat=B_all$GCV_mat
    TrHat=B_all$HatMatTrace
  }
  
  idxB=unlist(sapply(1:nList,function(ii)rep(ii,len=ncol(phiL[[ii]]))))
  
  if(method!='gq'){
    BL=split(B,idxB)
    out=unlist(mapply('%*%',phiL,BL,SIMPLIFY = F))  
  }
  
  if(method=='gq'){
    BL=split(B,idxB)
    out=mapply('%*%',phiL,BL,SIMPLIFY = F)
    outBL=lapply(ndimL,function(i)rep(0,i))
    for(iL in idxFD){
      outBL[[iL]][idxL[[iL]]]=out[[iL]]
    }
    outBL[-idxFD]=BL[-idxFD]
    out=unlist(outBL)
    xL=xL0
  }
  
  if(type=='cor'){
    corr=cor(yVec,drop(do.call(cbind,xL)%*%out))
    out=abs(corr)
    return(out)    
  }
  if(type=='a'){
    a=out/(sd(yVec)*abs(cor(yVec,do.call(cbind,xL)%*%out)))
    out=a
    return(out)    
  }
  if(type=='dir'){
    return(out)    
  }
  if(type=='all'){
    corr=abs(cor(yVec,do.call(cbind,xL)%*%out))
    a=out/(sd(yVec)*corr)
    return(list('corr'=corr,'a'=a,'K'=B_K,'gq'=gqWeight,
                'phiL'=phiL,'S'=B_S,'lam1'=B_lam1,'lam2'=B_lam2,
                'GCV_mat'=GCV_mat,'TraceHat'=TrHat))
  }
}




flars=function(x,y,method=c('basis','gq','raw'),max_selection,cv=c('gcv'),normalize=c('trace','rank','norm','raw'),lasso=TRUE,check=1,select=TRUE,VarThreshold=0.1,SignThreshold=0.8,control=list()){
  
  x=lapply(x,scale)
  xMean=lapply(x,function(i)attr(i,"scaled:center"))
  xSD=lapply(x,function(i)attr(i,"scaled:scale"))
  x=lapply(x,function(i){
    da=i
    attr(da,"scaled:center")=NULL
    attr(da,"scaled:scale")=NULL
    da
  })
  yMean=mean(y)
  ySD=sd(y)
  y=as.numeric(scale(y))
  nx=length(x)
  nsample=nrow(x[[1]])
  max_iter=max_selection
  alpha0=p2_norm0=NULL
  r=as.matrix(y)
  cv=cv[1]
  
  beta0=rep(0,ncol(do.call('cbind',x)))
  beta=betaM=NULL
  gamma0=mapply(function(a,b)as.matrix(rep(0,each=b)),rep(0,length(x)),sapply(x,ncol),SIMPLIFY = F)
  betaGQ=NULL
  
  
  A=A0=NULL
  method.=method
  if(cv=='gcv'){
    GCV=T

    cor1=sapply(x,function(ixL) fccaGen(xL = ixL,yVec=r[,1],type='cor',method=method.,GCV=T,control=control)) # find maximum correlation
  }
  
  idx1=as.numeric(which.max(cor1))
  A0=c(A0,idx1);A=as.numeric(na.omit(A0[A0>0]))
  iter=1
  varSplit=NULL
  
  MaxIter=min(max_iter,length(x))
  p0=NULL
  dfIter=NULL
  corIter=NULL
  SignCheckF0=NULL
  Checked=F
  cat('iter  ')
  #### start iteration ####
  while(iter <=MaxIter){
   cat(iter,' ')
   GCV0=1
    
    gamma_all=fccaGen(x[A],r[,iter,drop=F],type='all',method=method.,GCV=GCV0,control=control)
    
    gamma=split(gamma_all$a,unlist(mapply(rep,1:length(x[A]),each=sapply(x[A],ncol))))
    
    gamma2=gamma0
    gamma2[A]=gamma
    gamma2=as.matrix(unlist(gamma2))
    
    S0=gamma_all$S
    S0_norm=1
    
    p2=do.call('cbind',x)%*%gamma2*S0_norm
    p2_norm=1/sd(p2)
    p2=p2*p2_norm
    
    if(cv=='gcv'){
      GCV=T
      cor2=sapply(x[-A],function(ixL) fccaGen(xL = ixL,yVec=r[,iter],type='all',method=method.,GCV=T,control=control),simplify = F) # find maximum correlation
    }
    S=lapply(cor2,function(i)i$S)
    
    
    SidxNull=sapply(S,is.null)
    if(sum(unlist(SidxNull))!=0)
      S[SidxNull]=lapply(x[-A][SidxNull],function(i)tcrossprod(i)/drop(crossprod(i)))
    
    
    normalize=normalize[1]
    S_norm=sapply(S,function(i){
      if(normalize=='trace'){
        if(!is.null(i)) return(1/sum(diag(i))) #trace
      }
      if(normalize=='rank'){
        if(!is.null(i)) return(1/rankMatrix(i)) #rank
      }
      if(normalize=='norm'){
        if(!is.null(i)) return(1/(sqrt(sum(diag(crossprod(i)))))  ) #matrix norm
      }
      if(normalize=='raw'){
        if(!is.null(i)) return(1) #raw
      }
    })
    
    pBar=tcrossprod(p2)/drop(crossprod(p2))
    pk_norm=1
    if(normalize=='trace') pk_norm=1/sum(diag(pBar)) #trace
    
    if(normalize=='rank') pk_norm=1/rankMatrix(pBar) #rank
    
    if(normalize=='norm') pk_norm=1/(sqrt(sum(diag(cp(pBar))))) #matrix norm
    
    if(normalize=='raw') pk_norm==1
    
    SBar=mapply('*',S,S_norm,SIMPLIFY = F)
    pBar=pBar*pk_norm
    
    MatAlpha=matrix(0,ncol=2,nrow=length(x))
    colnames(MatAlpha)=c('plus','minus')
    MatAlpha[A,]=1e5
    
    Plus=Minus=rep(NA,length(cor2))
    
    aa=bb=cc=vector('list',length(cor2))
    
    try(aa<-lapply(SBar,function(i){
      cp(p2,i-pBar)%*%p2
    }),T)
    try(bb<-lapply(SBar,function(i){
      -2*cp(r[,iter],i-pBar)%*%p2
    }),T)
    try(cc<-lapply(SBar,function(i){
      cp(r[,iter],i-pBar)%*%r[,iter]
    }),T)
    
    aa=unlist(aa)
    bb=unlist(bb)
    cc=unlist(cc)
    
    bb2ac=bb^2-4*aa*cc
    bb2ac[bb2ac<0]=NA
    try(Plus<-(-bb+sqrt(bb2ac))/(2*aa),T)
    try(Minus<-(-bb-sqrt(bb2ac))/(2*aa),T)
    
    MatAlpha[-A,1]=unlist(Plus)
    MatAlpha[-A,2]=unlist(Minus)
    MatAlpha[MatAlpha<0]=NA
    alpha0=c(alpha0,as.vector(MatAlpha)[which.min(as.vector(MatAlpha))])
    ols=0
    if(alpha0[iter]==1e5){
      idx2=NA
      m0=lm(r[,iter]~0+p2)
      alpha0[iter]=m0$coefficients
      ols=1
    }
    if(alpha0[iter]!=1e5 & ols==0){
      idx2=c(which(MatAlpha[,1]==alpha0[iter]),which(MatAlpha[,2]==alpha0[iter]))  
    }

    A0=c(A0,idx2)
    A=as.numeric(na.omit(A0[A0>0]))
    
    beta=cbind(beta,as.vector(gamma2)*alpha0[iter]*p2_norm)
    p2_norm0=c(p2_norm0,p2_norm)
    
    p0=cbind(p0,p2)
    
    ## finish main calculation
    if(rmse(r[,iter],do.call(cbind,x)%*%beta[,iter,drop=F])>rmse(r[,iter],-do.call(cbind,x)%*%beta[,iter,drop=F])){
      alpha0[iter]=alpha0[iter]
      beta[,iter]=-beta[,iter]
    }
    
    if(iter==1)
      betaM=cbind(betaM,beta[,iter])
    if(iter>1)# & !Checked)
      betaM=cbind(betaM,betaM[,iter-1]+beta[,iter])
    # if(iter>1 & Checked)
    #   betaM=cbind(betaM,beta[,iter])
    
    ## start checking
    Checked==F
    if(lasso==T){
      if(check==1){
        betaL=split(betaM[,iter],unlist(mapply(rep,seq_along(x),each=sapply(x,ncol))))
        xSplit=mapply('%*%',x,betaL)
        tmpVar=apply(xSplit,2,var)
        varSplit=cbind(varSplit,tmpVar)
        VarIdx=as.logical(apply(varSplit,1,function(i){
          tail(i,1)<max(i[max(1,tail(which(i==0),1)):(length(i)-1)]) & tail(i,n = 1)>0 & tail(i,n = 1)<VarThreshold
        }))
        if(length(which(VarIdx==TRUE))>0){
          Checked=T
          RmIdx=which.min(tmpVar[VarIdx])
          RmIdx=which(VarIdx==TRUE)[RmIdx]
          RmIdxBeta=sum(sapply(betaL[1:RmIdx],length))-rev(seq_along(betaL[[RmIdx]]))+1
          betaM[RmIdxBeta,iter]=0
          A0[which(A0==RmIdx)]=A0[which(A0==RmIdx)]*(-1)-iter*0.0001
          A=as.numeric(na.omit(A0[A0>0]))
          
        }
      }
      
      
      
      if(check==2){
        betaL=split(betaM[,iter],unlist(mapply(rep,seq_along(x),each=sapply(x,ncol))))
        if(iter>2){
          SignCheckF=mapply(crossprod,lapply(betaL,sign),lapply(betaL0,sign))/sapply(betaL,function(iBetaL) length(iBetaL[iBetaL!=0]))
          SignCheckF0=cbind(SignCheckF0,as.numeric(SignCheckF))
          SignCheckIdx=which(SignCheckF<=SignThreshold & SignCheckF!=0)
          if(length(SignCheckIdx)>0){
            SignCheckIdx=SignCheckIdx[which.min(SignCheckF[SignCheckIdx])]
            A0[which(A0==SignCheckIdx)]=A0[which(A0==SignCheckIdx)]*(-1)-iter*0.0001
            A=as.numeric(na.omit(A0[A0>0]))
            BetaIdxL=split(1:nrow(betaM),unlist(mapply(rep,seq_along(x),each=sapply(x,ncol))))
            RmIdxBeta=unlist(BetaIdxL[SignCheckIdx])
            betaM[RmIdxBeta,iter]=0
          }
        }
        betaL0=betaL
      }
      
      ## finish checking
      
    }
    
    r=cbind(r,y-do.call(cbind,x)%*%betaM[,iter,drop=F])
    if(Checked & !select){
      r[,ncol(r)]=r[,1]
    }
    if(is.null(gamma_all$S)){
      gamma_all$S=tcrossprod(do.call(cbind,x[A]))
    }
    dfIter=c(dfIter,list(diag(1,nsample)-gamma_all$S))
    iter=iter+1
    corIter=c(corIter,cor(r[,iter],p2))
    p3=p2
    
  }
  betaM0=betaM
  
  Mu=Beta=NULL
  for(colI in 1:ncol(betaM0)){
    bb=betaM0[,colI,drop=F]
    b1=unlist(xMean)
    b2=bb/unlist(xSD)
    y1=yMean
    y2=ySD
    Mu=c(Mu,y1-y2*sum(b1*b2))
    Beta=cbind(Beta,y2*b2)
  }
  
  betaOLSdir=fccaGen(x,y,type='all')
  Sigma2Bar=var(y-do.call(cbind,x)%*%betaOLSdir$a)
  dfIterS=vector('list',length(dfIter))
  dfIterS[[1]]=dfIter[[1]]
  for(iDF in 2:length(dfIter)){
    dfIterS[[iDF]]=dfIter[[iDF]]%*%dfIterS[[iDF-1]]
  }
  dfIterS=lapply(dfIterS,function(iDFS) diag(1,nsample)-iDFS)
  dfIterS=sapply(dfIterS,function(iDFS) sum(diag(iDFS)))
  
  
  RSS=apply(do.call(cbind,x)%*%betaM0,2,function(cpI) sqrt(crossprod(y-cpI)))
  dfIterS=dfIterS[seq_along(RSS)]
  
  R2=1-RSS^2/crossprod(scale(y))
  R2Adj=1-(1-R2)*(nsample-1)/(nsample-dfIterS-1)
  
  
  Cp=RSS/Sigma2Bar-nsample+2*dfIterS
  CpTrue=NULL
  
  aicTrue=NULL
  aic=2*RSS/(2*Sigma2Bar)+2*nsample*(log(sqrt(Sigma2Bar))+log(2*pi))+dfIterS
  
  bicTrue=NULL
  bic=2*RSS/(2*Sigma2Bar)+2*nsample*(log(sqrt(Sigma2Bar))+log(2*pi))+dfIterS*log(nsample)
  
  StopStat=list(RSS=RSS,R2=R2,R2Adj=R2Adj,Cp=Cp,aic=aic,bic=bic,corIter=corIter)
  
  CD=abs(StopStat$corIter[-length(StopStat$corIter)]*alpha0[-1])
  
  
  out=list(Mu=Mu,
           Beta=Beta,
           # actual estimations
           alpha=alpha0,
           p2_norm=p2_norm0[1:(iter-1)],
           AllIndex=A0,
           # all the index
           index=A,
           # all the index except NA
           CD=CD,
           method=method,
           resid=r,
           RowMeans=xMean,
           RowSds=xSD,
           yMean=yMean,
           ySD=ySD,
           p0=p0,
           cor1=cor1,
           lasso=lasso,
           df=dfIterS,
           Sigma2Bar=Sigma2Bar,
           StopStat=StopStat,
           varSplit=varSplit,
           SignCheckF=SignCheckF0
  )
  class(out)='flars'
  return(out)
}


predict.flars=function(object,newdata,...){
  if (!inherits(object, "flars")) stop("Not a legitimate 'flars' object")
  
  xL=newdata
  if(!is.list(xL)) xL=list(xL)
  if(is.list(xL)) xL=lapply(xL,as.matrix)
  
  nList=length(xL)
  ndimL=sapply(xL,ncol)
  idxFD=seq_along(xL)[ndimL>1]
  
  if(length(idxFD)==0){
    x=sapply(xL,as.matrix)
    mu=object$Mu
    Beta=object$Beta
    
    yhattmp=matrix(rep(mu,each=nrow(x)),nrow=nrow(x),byrow = F)+x%*%Beta
  }
  
  if(length(idxFD)>0){
    mu=object$Mu
    Beta=object$Beta
    yhattmp=matrix(rep(mu,each=nrow(xL[[1]])),nrow=nrow(xL[[1]]),byrow = F)+do.call(cbind,xL)%*%Beta
  }
  
  return(yhattmp)
}


flars_TrainTest=function(seed=1,nsamples=120,nTrain=80,var_type=c('f','m'),VarThreshold0=0.1,SignThreshold0=0.8,cor_type=1:5,lasso=TRUE,check=1,uncorr=T,nVar=8,Discrete_Norm_ID=1:12,NoRaw_max=12,raw_max=9,hyper=NULL,RealX=NULL,RealY=NULL,dataL=NULL,nCor=0,control=list()){

  selection=Discrete_Norm_ID
  var_type=var_type[1]
  cor_type=cor_type[1]
  
  if(is.null(dataL) & is.null(RealX) & is.null(RealY)){
    dataL=data_generation(seed = seed,nsamples = nsamples,hyper = hyper,var_type =var_type,cor_type = cor_type,uncorr = uncorr,nVar = nVar)
  }
  
  fsm=dataL$x
  y=dataL$y 
  TrainIdx=seq(nTrain)
  TestIdx=seq(nsamples)[-TrainIdx]
  
  fsmTrain=lapply(fsm,function(fsmI) fsmI[TrainIdx,,drop=F])
  fsmTest=lapply(fsm,function(fsmI) fsmI[TestIdx,,drop=F])
  yTrain=y[TrainIdx]
  yTest=y[TestIdx]
  
  noise=dataL$noise

  if(!is.null(RealX) & !is.null(RealY)){
    fsm=RealX
    y=RealY
  }
  
  
  
  out=vector('list',length = 12)
  
  if(nCor<=1){
    iSam='seed'
    if(1%in%selection){
      test_Basis_Norm=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[1],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[2],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Basis_Norm
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Basis_Norm$yhat=cbind(yTest,yhattmp)
      
      out[[1]]=test_Basis_Norm
      cat(iSam,'.1',' ',sep = '')
    }
    
    if(2%in%selection){
      test_Basis_Rank=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[1],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[4],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Basis_Rank
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Basis_Rank$yhat=cbind(yTest,yhattmp)
      
      out[[2]]=test_Basis_Rank
      cat(iSam,'.2',' ',sep = '')
    }
    
    if(3%in%selection){
      test_Basis_Trace=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[1],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[1],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Basis_Trace
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Basis_Trace$yhat=cbind(yTest,yhattmp)
      
      out[[3]]=test_Basis_Trace
      cat(iSam,'.3',' ',sep = '')
    }
    
    if(4%in%selection){
      test_Basis_Raw=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[1],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[3],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Basis_Raw
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Basis_Raw$yhat=cbind(yTest,yhattmp)
      
      out[[4]]=test_Basis_Raw
      cat(iSam,'.4',' ',sep = '')
    }
    
    if(5%in%selection){
      test_GQ_Norm=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[2],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[2],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_GQ_Norm
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_GQ_Norm$yhat=cbind(yTest,yhattmp)
      
      out[[5]]=test_GQ_Norm
      cat(iSam,'.5',' ',sep = '')
    }
    
    if(6%in%selection){
      test_GQ_Rank=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[2],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[4],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_GQ_Rank
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_GQ_Rank$yhat=cbind(yTest,yhattmp)
      
      out[[6]]=test_GQ_Rank
      cat(iSam,'.6',' ',sep = '')
    }
    
    if(7%in%selection){
      test_GQ_Trace=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[2],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[1],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_GQ_Trace
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_GQ_Trace$yhat=cbind(yTest,yhattmp)
      
      out[[7]]=test_GQ_Trace
      cat(iSam,'.7',' ',sep = '')
    }
    
    if(8%in%selection){
      test_GQ_Raw=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[2],max_selection = NoRaw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[3],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_GQ_Raw
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_GQ_Raw$yhat=cbind(yTest,yhattmp)
      
      out[[8]]=test_GQ_Raw
      cat(iSam,'.8',' ',sep = '')
    }
    
    if(9%in%selection){
      test_Raw_Norm=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[3],max_selection = raw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[2],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Raw_Norm
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Raw_Norm$yhat=cbind(yTest,yhattmp)
      
      out[[9]]=test_Raw_Norm
      cat(iSam,'.9',' ',sep = '')
    }
    
    if(10%in%selection){
      test_Raw_Rank=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[3],max_selection = raw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[4],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Raw_Rank
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Raw_Rank$yhat=cbind(yTest,yhattmp)
      
      out[[10]]=test_Raw_Rank
      cat(iSam,'.A',' ',sep = '')
    }
    
    if(11%in%selection){
      test_Raw_Trace=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[3],max_selection = raw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[1],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Raw_Trace
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Raw_Trace$yhat=cbind(yTest,yhattmp)
      
      out[[11]]=test_Raw_Trace
      cat(iSam,'.B',' ',sep = '')
    }
    
    if(12%in%selection){
      test_Raw_Raw=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[3],max_selection = raw_max,lasso=lasso,normalize = c('trace','norm','raw','rank')[3],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      test=test_Raw_Raw
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test_Raw_Raw$yhat=cbind(yTest,yhattmp)
      
      out[[12]]=test_Raw_Raw
      cat(iSam,'.C',' ',sep = '')
    }
  }
  
  
  
  if(nCor>1){
    SelectionIdx=cbind(rep(1:3,each=4),rep(1:4,3))
    CatIdx=c(1:9,'A','B','C')
    res=mclapply(selection,function(iSam){
      methodIdx=SelectionIdx[iSam,1]
      normIdx=SelectionIdx[iSam,2]
      if(methodIdx<3) No_Max=NoRaw_max
      if(methodIdx==3) No_Max=raw_max
      test=flars(fsmTrain,yTrain,method = c('basis','gq','raw')[methodIdx],max_selection = No_Max,lasso=lasso,normalize = c('trace','norm','raw','rank')[normIdx],check=check,VarThreshold=VarThreshold0,SignThreshold=SignThreshold0,control=control)
      mu=test$Mu
      Beta=test$Beta
      yhattmp=matrix(rep(mu,each=length(TestIdx)),nrow=length(TestIdx),byrow = F)+do.call(cbind,fsmTest)%*%Beta
      test$yhat=cbind(yTest,yhattmp)
      
      cat('seed.',CatIdx[iSam],sep = '')
      return(test)
    },mc.cores = nCor,mc.set.seed=F)
    out[selection]=res
  }
  names(dataL$x)=paste0('x',seq_along(dataL$x))
  out$dataL=dataL
  names(out)=c('out_Basis_Norm','out_Basis_Rank','out_Basis_Trace','out_Basis_Raw','out_GQ_Norm','out_GQ_Rank','out_GQ_Trace','out_GQ_Raw','out_Raw_Norm','out_Raw_Rank','out_Raw_Trace','out_Raw_Raw','dataL')
  return(out)
}


fccaXXcv=function(xL1,xL2,method=c('basis','gq','raw'),centre = TRUE,tol=1e-7,Control1=list(),Control2=list(),alpha=10^seq(-6,1,len=10)){
  out=rep(0,length(alpha))
  for(i in seq_along(alpha)){
    #cat(i)
    Control1$pen1=c(Control1$pen1,alpha[i])[1]
    Control2$pen1=c(Control2$pen1,alpha[i])[1]
    Control2$pen2=Control2$pen2=0
    cvL1=cvL2=rep(0,nrow(xL1[[1]]))
    for(j in 1:nrow(xL1[[1]])){
      cvL1[j]= cvL2[j] =0
      xL1nI=lapply(xL1,function(ii) ii[-j,])
      xL2nI=lapply(xL2,function(ii) ii[-j,])
      xL1I=lapply(xL1,function(ii) ii[j,,drop=F])
      xL2I=lapply(xL2,function(ii) ii[j,,drop=F])
      
      tmp=fccaXX(xL1nI,xL2nI,method=method,centre = centre,control1=Control1,tol = tol,control2=Control2)
      tmp
      try(cvL1[j]<-Reduce('+',mapply('%*%',xL1I,lapply(tmp$coef1,function(ii) ii[,1,drop=F] ))),T)
      try(cvL2[j]<-Reduce('+',mapply('%*%',xL2I,lapply(tmp$coef2,function(ii) ii[,1,drop=F] ))),T)
      if(sum(c(cvL1[j], cvL2[j])%in%0)>0)
        break
    }
    if(length(cvL1)!=nrow(xL1[[1]]))
      out[i]=NA
    if(length(cvL1)==nrow(xL1[[1]]))
      out[i]=suppressWarnings(cor(cvL1,cvL2))
    Control1$pen1=Control2$pen1=NULL
  }
  idx_na=which(is.na(out))
  outcome=list(cor=out[-idx_na],alpha=alpha[-idx_na])
  return(outcome)
}



fccaPre=function(xL,centre=TRUE,method=c('basis','gq','raw'),control=list()){
  if(!is.list(xL)) xL=list(xL)
  
  if(is.list(xL)) xL=lapply(xL,as.matrix)
  
  nList=length(xL)
  ndimL=sapply(xL,ncol)
  idxFD=seq_along(xL)[ndimL>1]
  
  if(centre){
    xL=lapply(xL,scale)
    xMean=lapply(xL,function(i)attr(i,"scaled:center"))
    xSD=lapply(xL,function(i)attr(i,"scaled:scale"))
    xL=lapply(xL,function(i){
      da=i
      attr(da,"scaled:center")=NULL
      attr(da,"scaled:scale")=NULL
      da
    })
  }
  
  
  method=method[1]
  gqWeight=NULL
  
  
  pL=phiL=vector('list',length=nList)
  
  if(method=='basis'){
    con=list(nbasis=35,norder=6,pen1=10^(seq((-20),5,len=41)),pen2=0.01,t=seq(0,1,len=max(sapply(xL,ncol),na.rm = T)))
    con[(namc <- names(control))] <- control
    
    nbasis=con$nbasis;norder=con$norder;pen1=con$pen1
    pen2=con$pen2;t=con$t
    
    if(!is.list(t)){
      tRange=range(t)
      t=lapply(xL,function(i){
        if(ncol(i)==length(t)){
          return(t)
        }
        if(ncol(i)!=length(t) & ncol(i)>1){
          return(seq(tRange[1],tRange[2],len=ncol(i)))
        }
        if(ncol(i)==1){
          return((tRange[1]+tRange[2])/2)
        }
      })
    }
    if(is.list(t)){
      if(length(t)!=nList)
        t=lapply(xL,function(i)return(t[[1]]))
    }
    
    
    if(length(nbasis)!=length(norder)){
      nbasis=nbasis[1]
      norder=norder[1]
      warning('length of nbasis is different from length of norder')
    }
    if(length(nbasis)!=nList)
      nbasis=rep(nbasis[1],nList)
    if(length(norder)!=nList)
      norder=rep(norder[1],nList)
    
    for(iL in idxFD){
      tRange=range(t[[iL]])
      nb=nbasis[iL];no=norder[iL]
      tiL=t[[iL]]
      spline=create.bspline.basis(tRange,nbasis = nb,norder = no)
      MS=eval.basis(tiL,spline,0)
      MS2=eval.basis(tiL,spline,2)
      
      pL[[iL]]=crossprod(MS2)
      phiL[[iL]]=MS/(ncol(xL[[iL]]))
    }
  }
  
  if(method=='gq'){
    con=list(nP=18,pen1=10^(seq((-20),5,len=21)),pen2=0.01,t=seq(-1,1,len=max(sapply(xL,ncol),na.rm = T)))
    con[(namc <- names(control))] <- control
    
    t=con$t;nP=con$nP;pen1=con$pen1;pen2=con$pen2
    
    if(!is.list(t)){
      tRange=range(t)
      t=lapply(xL,function(i){
        if(ncol(i)==length(t)){
          return(t)
        }
        if(ncol(i)!=length(t) & ncol(i)>1){
          return(seq(tRange[1],tRange[2],len=ncol(i)))
        }
        if(ncol(i)==1){
          return((tRange[1]+tRange[2])/2)
        }
      })
    }
    if(is.list(t)){
      if(length(t)!=nList)
        t=lapply(xL,function(i)return(t[[1]]))
    }
    if(GQ==0){
      .GlobalEnv$modelChronic  <-  readRDS(
        file.path(system.file('inst/extdata',package='flars'), 'GQ.rds'))
    }
    gqWeight=data.matrix(GQ[GQ$n==nP,c(1,2)])
    gqWeight2=data.matrix(GQ[GQ$n==(nP-2),c(1,2)])
    xi=as.numeric(gqWeight[,1])
    wi=as.matrix(gqWeight[,2])
    wi2=as.matrix(gqWeight2[,2])
    for(iL in idxFD){
      dimFx=dim(xL[[iL]])
      ncols=dimFx[2]
      nrows=dimFx[1]
      zi=t[[iL]]
      idx=as.numeric(sapply(xi,function(i) which.min(abs(zi-i))))
      xGQ=t(apply(xL[[iL]][,sort(idx)],1,function(i)i*(wi[order(idx)])))
      xL[[iL]]=xGQ
      
      pL[[iL]]=t(get2L(ncol(xL[[iL]]),diff(sort(xi))))%*%diag(wi2[,1])%*%get2L(ncol(xL[[iL]]),diff(sort(xi)))
      phiL[[iL]]=diag(wi[,1])#/(ncol(xL[[iL]])^2)
    }
  }
  
  if(method=='raw'){
    con=list(pen1=10^(seq((-10),-5,len=21)),pen2=0.01,t=seq(-1,1,len=max(sapply(xL,ncol),na.rm = T)))
    con[(namc <- names(control))] <- control
    t=con$t;pen1=con$pen1;pen2=con$pen2
    
    if(!is.list(t)){
      tRange=range(t)
      t=lapply(xL,function(i){
        if(ncol(i)==length(t)){
          return(t)
        }
        if(ncol(i)!=length(t) & ncol(i)>1){
          return(seq(tRange[1],tRange[2],len=ncol(i)))
        }
        if(ncol(i)==1){
          return((tRange[1]+tRange[2])/2)
        }
      })
    }
    if(is.list(t)){
      if(length(t)!=nList)
        t=lapply(xL,function(i)return(t[[1]]))
    }
    
    for(iL in idxFD){
      pL[[iL]]=crossprod(get2L(ncol(xL[[iL]]),diff(t[[iL]])))
      phiL[[iL]]=diag(1,ncol(xL[[iL]]))/(ncol(xL[[iL]]))
    }
  }
  
  if(length(idxFD)<nList){
    pL[-idxFD]=lapply(pL[-idxFD],function(i) i=matrix(0))
    phiL[-idxFD]=lapply(phiL[-idxFD],function(i) i=matrix(1))
  }
  W=phiL;W2=pL
  
  v=do.call(cbind,mapply('%*%',xL,W,SIMPLIFY = F))
  vv=crossprod(v)
  Pen=matrix(0,ncol=ncol(vv),nrow=nrow(vv))
  
  DimPhL=sapply(phiL,ncol)
  DimP1=DimPhL[1]
  DimPhL=DimPhL-DimP1
  for(i in seq_along(phiL)){
    Pen[1:(DimPhL[i]+DimP1)+DimP1*(i-1)+sum(DimPhL[1:(i-1)]),
        1:(DimPhL[i]+DimP1)+DimP1*(i-1)+sum(DimPhL[1:(i-1)])]=pen1[1]*W2[[i]]+pen2[1]*crossprod(W[[i]])
  }
  P=vv+Pen
  output=list(P=P,Xv=v,W=W,W2=W2)
  return(output)
}

fccaXX=function(xL1,xL2,centre=TRUE,method=c('basis','gq','raw'),control1=list(),control2=list(),tol=1e-7){
  pre1=fccaPre(xL1,centre=centre,method = method[1],control = control1)
  pre2=fccaPre(xL2,centre=centre,method = method[1],control = control2)
  
  M1=ginv(pre1$P)%*%crossprod(pre1$Xv,pre2$Xv)%*%ginv(pre2$P)%*%crossprod(pre2$Xv,pre1$Xv)
  M2=ginv(pre2$P)%*%crossprod(pre2$Xv,pre1$Xv)%*%ginv(pre1$P)%*%crossprod(pre1$Xv,pre2$Xv)
  
  out1=eigen(M1)
  out2=eigen(M2)
  
  rho=rho1=Re(out1$values)
  rho2=Re(out2$values)
  if(rho1[1]-rho2[1]>1e-7) warning("eigenvalues don't match")
  idx=rho>tol
  rho=rho[idx]
  a=Re(out1$vectors)[,idx,drop=F]
  b=Re(out2$vectors)[,idx,drop=F]
  if(length(idx)==0){
    rho=0
    coef1=coef2=1
  }
  if(length(idx)!=0){
    aDim=c(0,cumsum(sapply(pre1$W,ncol)))
    coef1=vector('list',length = length(xL1))
    for(i in seq_along(xL1)){
      coef1[[i]]=pre1$W[[i]]%*%a[(1+aDim[i]):aDim[i+1],]
    }
    bDim=c(0,cumsum(sapply(pre2$W,ncol)))
    coef2=vector('list',length = length(xL2))
    for(i in seq_along(xL2)){
      coef2[[i]]=pre2$W[[i]]%*%b[(1+bDim[i]):bDim[i+1],]
    }
  }
  
  out=list(corr=rho,coef1=coef1,coef2=coef2,corr2=rho2)
  return(out)
}
