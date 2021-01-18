FASTmrMLM<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svrad,svmlod,Genformat,CLO){

inputform<-Genformat
svlod<-svmlod

if(is.null(kk)){
  
  if(is.null(gen)==TRUE)
  {
    warning("Please input correct genotype dataset !")
  }else{
    XX1<-t(gen)
    x<-XX1[3:nrow(XX1),]
    rownames(x)<-NULL
    colnames(x)<-NULL
    X1<-as.matrix(x)
    rm(x,XX1)
    gc()
    n<-nrow(X1)
    m<-ncol(X1)
    ########kinship##########
    #kk1<-(X1%*%t(X1))/m
    kk1<-mrMLM::multiplication_speed(X1,t(X1))/m
    kk<-as.matrix(kk1) 
  }
  rm(kk1,X1)
  gc()
} 

if(is.null(psmatrix)){
  flagps<-1
}else{
  flagps<-0
}

if(is.null(svpal)==TRUE||is.null(svrad)==TRUE||is.null(svlod)==TRUE){
  warning("Please set parameter!")
}

if((svpal<0)||(svpal>1))
{
  warning("Please input critical P-value between 0 and 1!")
}
if(svrad<0)
{
  warning("Please input search radius (kb) of candidate gene: > 0 !")
}
if(svlod<0)
{
  warning("Please input critical LOD score: > 0 !")
}

if(exists("gen")==FALSE)
{
  warning("Please input correct genotype dataset !")
}
if(exists("phe")==FALSE)
{
  warning("Please input correct phenotype dataset !")
}
if(exists("kk")==FALSE)
{
  warning("Please input correct kinship (K) dataset !")
}
if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
{
  warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset!")
}

if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2)))&&(svpal>=0)&&(svpal<=1)&&(svrad>0)&&(svmlod>=0))
{

parmsShow<-NULL
wan<-NULL
parms<-NULL
parms.pchange<-NULL

multinormal<-function(y,mean,sigma)
{
  pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
  return (pdf_value)
}

ebayes_EM<-function(x,z,y)
{
  n<-nrow(z);k<-ncol(z)
  
  if(abs(min(eigen(crossprod(x,x))$values))<1e-6){
    b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
  }else{
    b<-solve(crossprod(x,x))%*%(crossprod(x,y))
  }
  
  v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
  u<-matrix(rep(0,k),k,1)
  v<-matrix(rep(0,k),k,1)
  s<-matrix(rep(0,k),k,1)
  for(i in 1:k)
  {
    zz<-z[,i]
    s[i]<-((crossprod(zz,zz)+1e-100)^(-1))*v0
    u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
    v[i]<-u[i]^2+s[i]
  }
  
  vv<-matrix(rep(0,n*n),n,n);
  for(i in 1:k)
  {
    zz<-z[,i]
    vv=vv+tcrossprod(zz,zz)*v[i]
  }
  vv<-vv+diag(n)*v0
  
  iter<-0;err<-1000;iter_max<-500;err_max<-1e-8
  tau<-0;omega<-0
  while((iter<iter_max)&&(err>err_max))
  {
    iter<-iter+1
    v01<-v0
    v1<-v
    b1<-b
    vi<-solve(vv)
    xtv<-crossprod(x,vi)
    
    if(ncol(x)==1)
    {
      b<-((xtv%*%x)^(-1))*(xtv%*%y)
    }else{
      if(abs(min(eigen(xtv%*%x)$values))<1e-6){
        b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
      }else{
        b<-solve(xtv%*%x)%*%(xtv%*%y)
      }
    }
    r<-y-x%*%b
    ss<-matrix(rep(0,n),n,1)
    for(i in 1:k)
    {
      zz<-z[,i]
      zztvi<-crossprod(zz,vi)
      u[i]<-v[i]*zztvi%*%r
      s[i]<-v[i]*(1-zztvi%*%zz*v[i])
      v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
      ss<-ss+zz*u[i]
    }
    v0<-as.numeric(crossprod(r,(r-ss))/n)
    
    vv<-matrix(rep(0,n*n),n,n)
    for(i in 1:k)
    {
      zz<-z[,i]
      vv<-vv+tcrossprod(zz,zz)*v[i]
    }
    vv<-vv+diag(n)*v0
    
    err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
    beta<-t(b)
    sigma2<-v0
  }
  
  wang<-matrix(rep(0,k),k,1)
  for (i in 1:k){
    stderr<-sqrt(s[i]+1e-20)
    t<-abs(u[i])/stderr
    f<-t*t
    p<-pchisq(f,1,lower.tail = F)
    wang[i]<-p
  }
  
  return(list(u=u,sigma2=sigma2,wang=wang))
}

likelihood<-function(xxn,xxx,yn,bbo)
{
  nq<-ncol(xxx)
  ns<-nrow(yn)
  at1<-0
  
  if(is.null(bbo)==TRUE){
    ww1<-1:ncol(xxx)
    ww1<-as.matrix(ww1)
  }else{
    ww1<-as.matrix(which(abs(bbo)>1e-5))
  }
  at1<-dim(ww1)[1]
  lod<-matrix(rep(0,nq),nq,1)
  if(at1>0.5)
    ad<-cbind(xxn,xxx[,ww1])
  else
    ad<-xxn
  if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
    bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
  else
    bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
  vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
  ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
  sub<-1:ncol(ad);
  if(at1>0.5)
  {
    for(i in 1:at1)
    {
      ij<-which(sub!=sub[i+ncol(xxn)])
      ad1<-ad[,ij]
      if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
        bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
      else
        bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn) 
      vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
      ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
      lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
    }
  }
  return (lod)
}

mixed1<-function(xu,yu,theta1){
  
  loglike<-function(theta1){
    lambda<-exp(theta1)
    logdt<-sum(log(lambda*delta+1))
    h<-1/(lambda*delta+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,q,1)
    xx<-matrix(0,q,q)
    for(i in 1:q){
      yx[i]<-sum(yu*h*xu[,i])
      for(j in 1:q){
        xx[i,j]<-sum(xu[,i]*h*xu[,j])
      }
    }
    loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))-0.5*(n-q)
    return(-loglike)  
  }
  grad<-function(theta1){
    lambda<-exp(theta1)
    h<-1/(lambda*delta+1)
    d<-diag(delta,nrow(X1),nrow(X1))
    hinv<-diag(1/(lambda*delta+1),nrow(X1),nrow(X1))
    yy<-sum(yu*h*yu)
    yx<-matrix(0,q,1)
    xx<-matrix(0,q,q)
    for(i in 1:q){
      yx[i]<-sum(yu*h*xu[,i])
      for(j in 1:q){
        xx[i,j]<-sum(xu[,i]*h*xu[,j])
      }
    }
    pp=hinv- hinv%*%xu%*%solve(xx)%*%t(xu)%*%hinv
    sigma<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)		
    f= -0.5*{sum(diag(pp%*%d))-1/sigma*(t(yu)%*%pp%*%d%*%pp%*%yu)}
    return(c(-f))
  }
  parm<-optim(par=theta,fn=loglike,gr=grad,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
  lambda<-(parm$par)
  return(c(lambda))
}

lll<- function(theta){
  lambdak<-exp(theta)
  deth<-1+lambdak*g1
  tmp<-lambdak*1/deth
  yHy<-yy-zy%*%tmp%*%zy
  yHx<-yx-zx%*%tmp%*%zy
  xHx<-xx-zx%*%tmp%*%t(zx)
  logdt2<-log(deth)
  ll<- -0.5*logdt2-0.5*(n-q)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
  return(-ll)
}
grad2<- function(theta){
  lambdak<-exp(theta)
  deth<-1+lambdak*g1
  tmp<-lambdak*1/deth
  yHy<-yy-zy%*%tmp%*%zy
  yHx<-yx-zx%*%tmp%*%zy
  xHx<-xx-zx%*%tmp%*%t(zx)
  zHy<-zy-zz%*%tmp%*%zy
  zHx<-zx-zx%*%tmp%*%zz
  zHz<-zz-zz%*%tmp%*%zz
  sigma2<-(yHy-t(yHx)%*%solve(xHx)%*%yHx)/(n-q)
  f<- -0.5*{(zHz-t(zHx)%*%solve(xHx)%*%zHx)-(zHy-t(zHx)%*%solve(xHx)%*%yHx)^2/sigma2}
  return(c(-f))
}


fixed2<-function(lambdak){
  deth<-1+lambdak*g1
  tmp<-lambdak*1/deth
  yHy<-yy-zy%*%tmp%*%zy
  yHx<-yx-zx%*%tmp%*%zy
  xHx<-xx-zx%*%tmp%*%t(zx)
  zHy<-zy-zz%*%tmp%*%zy
  zHx<-zx-zx%*%tmp%*%zz
  zHz<-zz-zz%*%tmp%*%zz
  beta<-solve(xHx,yHx)
  tmp2<-solve(xHx)
  sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-q)
  gamma<-lambdak*zHy-lambdak*t(zHx)%*%tmp2%*%yHx
  var<-abs((lambdak*diag(1)-lambdak*zHz*lambdak)*as.numeric(sigma2))
  wald<-gamma^2/var
  stderr<-sqrt(diag(var))
  p_value<-pchisq(wald,1,lower.tail = F) 
  result<-list(gamma,stderr,beta,sigma2,p_value,wald)
  return(result)
}

y<-as.matrix(phe)
XX1<-t(gen)

x<-XX1[3:nrow(XX1),]
rownames(x)<-NULL
colnames(x)<-NULL

X1<-as.matrix(x)
rm(x)
gc()
n<-nrow(X1)
m<-ncol(X1)
########kinship##########
xxx<-matrix(1,n,1)
xxx<-matrix()
if (is.null(psmatrix)==TRUE)
{
  xxx<-matrix(1,n,1)
}else{
  ps<-as.matrix(psmatrix)
  xxx<-cbind(matrix(1,n,1),ps) 
}

qq<-eigen(kk)
delta<-qq[[1]]
d<-diag(delta,n,n)
uu<-qq[[2]]
q<-ncol(xxx)
waving<-svrad
xu<-t(uu)%*%xxx
zkk<-t(uu)%*%X1
theta1<-0
theta<-0 

rm(kk,d,qq)
gc()

ll<-numeric()
y<-as.matrix(y)
yu<-t(uu)%*%y
ll<-numeric()
omeg<-mixed1(xu,yu,theta1)
delta1<-1/sqrt(delta*exp(omeg)+1)
d1<-diag(delta1,nrow(X1),nrow(X1))
yc<-d1%*%yu
yy<-sum(yc*1*yc)
xc<-d1%*%xu
yx<-matrix(0,q,1)
for(i in 1:q){
  yx[i]<-sum(yc*1*xc[,i])
}
binv<-diag(1,nrow(X1),nrow(X1))
xx<-matrix(0,q,q)
for(i in 1:q){
  for(j in 1:q){
    xx[i,j]<-sum(xc[,i]*1*xc[,j])
  }
}
zkk1<-d1%*%zkk

rm(d1,uu)
gc()

cl.cores <- detectCores()
if((cl.cores<=2)||(is.null(CLO)==FALSE)){
  cl.cores<-1
}else if(cl.cores>2){
  if(cl.cores>10){
    cl.cores<-10
  }else {  
    cl.cores <- detectCores()-1
  }
}   
cl <- makeCluster(cl.cores)
registerDoParallel(cl)  
mat=foreach(j=1:m, .multicombine=TRUE, .combine = 'rbind')%dopar%
{
  
  zc<-as.matrix(zkk1[,j])
  uu1<-as.matrix(zc)%*%t(as.matrix(zc))
  g1<-sum(diag(uu1))
  zy<-as.matrix(sum(yc*1*zc))
  zz<-as.matrix(sum(zc*1*zc))
  zx<-matrix(0,q,1)
  for(i in 1:q){
    zx[i]<-sum(xc[,i]*1*zc)
  }
  par<-tryCatch(optim(par=theta,fn=lll,hessian = TRUE,gr=grad2,method="L-BFGS-B",lower=-10,upper=10), error=function(e) optim(par=theta,fn=lll,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10))
  lambda<-exp(par$par)
  conv<-par$convergence
  fn1<-par$value
  hess<-par$hessian
  parmfix<-fixed2(lambda)
  gamma<-parmfix[[1]]
  stderr<-parmfix[[2]]
  beta<-parmfix[[3]][1,]
  sigma2<-parmfix[[4]]
  p_wald<-parmfix[[5]]
  sigma2g<-lambda*sigma2
  wald<-parmfix[[6]]
  fn0<-lll(c(-Inf))
  lrt<-2*abs(fn0-fn1)
  p_lrt<-pchisq(lrt,1,lower.tail = F)
  parm0<-c(j,beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)
}
stopCluster(cl)

rm(zkk,zkk1)
gc()

ll<-rbind(ll,mat)
parms1<-as.matrix(ll)
rownames(parms1)<-NULL
newparm<-cbind(gen[,1:2],parms1[,2:8])
parms<-newparm
parms.pchange<-parms
parmsp<-as.matrix(parms.pchange[,9])
locsub<-which(parmsp==0)
if(length(locsub)!=0){
  pmin<-min(parmsp[parmsp!=0])
  subvalue<-10^(1.1*log10(pmin))
  parms.pchange[locsub,9]<-subvalue
}else{
  parms.pchange<-parms
}

if(inputform==1){
  #output result1 using mrMLM numeric format
  parmsShow<-parms
  tempparms<-parms[,3:9]
  tempparms[,7]<--log10(tempparms[,7])
  tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
  tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
  parmsShow<-cbind(genRaw[-1,1],parms[,1:2],tempparms,genRaw[-1,4])
  colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (FASTmrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (FASTmrMLM)'","Genotype for code 1")
  
}
if(inputform==2){
  #output result1 using mrMLM character format
  parmsShow<-parms
  outATCG<-matrix(outATCG,,1)
  tempparms<-parms[,3:9]
  tempparms[,7]<--log10(tempparms[,7])
  tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
  tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
  parmsShow<-cbind(genRaw[-1,1],parms[,1:2],tempparms,outATCG)
  colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (FASTmrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (FASTmrMLM)'","Genotype for code 1")
  
}
if(inputform==3){
  #output result1 using TASSEL format
  parmsShow<-parms
  outATCG<-matrix(outATCG,,1)
  outATCG<-unlist(strsplit(outATCG,""))
  outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
  tempparms<-parms[,3:9]
  tempparms[,7]<--log10(tempparms[,7])
  tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
  tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
  parmsShow<-cbind(genRaw[-1,1],parms[,1:2],tempparms,outATCG)
  colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","Mean","Sigma2","Sigma2_k","SNP effect (FASTmrMLM)","Sigma2_k_posteriori","Wald","'-log10(P) (FASTmrMLM)'","Genotype for code 1")
  
}



p<-as.vector(parms1[,8])
ans<-p.adjust(p, method = "bonferroni", n = length(p))

rm(gen)
gc()

##########p is parameter########
sigg<-as.vector(which(p<=svpal))
le1<-length(sigg)

if(le1!=0){
  
  if (length(which(ans<0.05))!=0)
  {
    siggbh<-which(ans<0.05)
    nnn1<-cbind(XX1[1,],XX1[2,])
    setloci<-siggbh
    setposi<-c(XX1[2,siggbh])
    num<-dim(nnn1)[1]
    endresult<-numeric()
    for (t in 1:length(siggbh))
    {
      for (i in 1:num){
        temp<-numeric()
        if ((XX1[1,i]==XX1[1,(setloci[t])])&&(abs(nnn1[i,2]-setposi[t])<=waving))
        {
          temp<-cbind(matrix(nnn1[i,],1,),i)
          endresult<-rbind(endresult,temp)
        }
      }
    }
    end<-as.vector(endresult[,3])
    sigg2<-sigg[!sigg%in% end]
    sigg1<-sort(c(siggbh,sigg2))
  }else{
    sigg1<-sigg
  }
  if (length(sigg1)>nrow(X1))
  {
    larsres<-lars(X1[,sigg1], y, type = "lar",trace = FALSE, normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps, use.Gram=FALSE) 
    larsc2<-sigg1[which(larsres$entry!=0)]
    if(length(which(larsres$entry>nrow(X1)))!=0)
    {
      ad1<-sigg1[which(larsres$entry>nrow(X1))]
      larsc<-larsc2[!larsc2%in%ad1]
    }else{
      larsc<-larsc2
    }
  }else{
    larsc<-sigg1 
  }
  
  z<-matrix(1,nrow(X1),1)
  z<-matrix()
  
  if (is.null(psmatrix)==TRUE)
  {
    z<-matrix(1,nrow(X1),1)
  }else{
    z<-cbind(matrix(1,nrow(X1),1),psmatrix) 
  }
  le1<-length(larsc)
  xxxnew11<-as.matrix(X1[,larsc])
  
  u1<-ebayes_EM(z,xxxnew11,y)
  obj<-u1$u
  result1<-matrix(0,m,1)
  for (i in 1:le1)
  {
    result1[(larsc)[i],1]=obj[i]
  } 
  Res<- t(as.matrix((rowSums(result1)/ncol(result1))))
  Res1<-as.vector(Res)	
  sig1<-which(abs(Res1)>=1e-5)
  le2<-length(which(abs(Res1)>=1e-5))
  
  if(le2!=0){
    bbo<-matrix(0,le2,1)
    for (i in 1:le2){
      bbo[i,]=Res1[sig1[i]]
    }
    xxxx<-as.matrix(X1[,sig1])
    yn<-as.matrix(y)
    xxn<-z
    lod<-likelihood(xxn,xxxx,yn,bbo)
    
    her1<-vector(length=le2)
    for (i in 1:le2){
      p1<-length(as.vector(which(X1[,sig1[i]]==1)))/length(X1[,sig1[i]])
      p2<-1-p1
      her1[i]=((p1+p2)-(p1-p2)^2)*(Res1[sig1[i]])^2
    }
    
    if(var(y)>=sum(her1)+u1$sigma2){
      her<-(her1/as.vector(var(y)))*100  
      
    }else{
      her<-(her1/(sum(her1)+u1$sigma2))*100 
    }
    
    slod<-cbind(sig1,lod,her)
    
    if(length(which(slod[,2]>=svlod))>=1){
      
      if(length(which(slod[,2]>=svlod))==1){
        sslod<-t(as.matrix(slod[which(slod[,2]>=svlod),]))
        sig1<-slod[which(slod[,2]>=svlod),1]
      }else if(length(which(slod[,2]>=svlod))>1){
        sslod<-slod[which(slod[,2]>=svlod),]
        sig1<-sslod[,1]
      }
      xxxx<-as.matrix(X1[,sig1])
      lod<-sslod[,2]  
      her<-sslod[,3]
      
      ii<-as.vector(sig1)
      qqq<-matrix(0,nrow=length(ii),ncol=6)
      qqq[,1]=as.matrix(ii)
      for (j in 1:length(ii)){
        qqq[j,2]=XX1[1,ii[j]]
        qqq[j,3]=XX1[2,ii[j]]
        qqq[j,4]=result1[ii[j],]
        
        qqq[j,5]=lod[j]
        qqq[j,6]=her[j]
      }
      
      rm(XX1,X1)
      gc()
      
      id<-which(qqq[,5]==0)
      
      if(length(id)!=dim(qqq)[1]){
        
        if(length(id)!=0){
          qqq1<-qqq[-id,]
        }else{
          qqq1<-qqq
        }
        xxmaf<-t(xxxx)
        leng.maf<-dim(xxmaf)[2]
        
        maf.fun<-function(snp){
          leng<-length(snp)
          snp1<-length(which(snp==1))
          snp11<-length(which(snp==-1))
          snp0<-length(which(snp==0))
          ma1<-(2*snp1+snp0)/(2*leng)
          ma2<-(2*snp11+snp0)/(2*leng)
          maf<-min(ma1,ma2)
          return(maf)
        }
        
        maf<-apply(xxmaf,1,maf.fun)
        maf<-as.matrix(round(maf,4))
        vee<-round(u1$sigma2,4)
        pee<-round(var(y),4)
        
        if(nrow(qqq1)>1){
          result<-as.matrix(qqq1[,-1])
          vees<-matrix("",nrow = nrow(result),1)
          pees<-matrix("",nrow = nrow(result),1)
          pees[1,1]<-pee
          vees[1,1]<-vee
          
        }else{
          result<-t(as.matrix(qqq1[,-1]))
          pees<-as.matrix(pee)
          vees<-as.matrix(vee)
        }
        
        
        if(nrow(qqq1)>1){
          result<-as.matrix(qqq1[,-1])
          result<-result
          temp<-as.matrix(result[,3:5])
          temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
          temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
          wan<-cbind(result[,1:2],temp)
          snp<-parmsShow[,11]
          
        }else{
          result<-t(as.matrix(qqq1[,-1]))
          result<-result
          temp<-t(as.matrix(result[,3:5]))
          temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
          temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
          wan<-cbind(t(as.matrix(result[,1:2])),temp)
          snp<-parmsShow[,11]
          
        }
        
        if(inputform==1){
          genRaw<-as.data.frame(genRaw)
          genraw<-genRaw[-1,1:4]
          
          wan_len<-dim(wan)[1]
          marker<-character()
          snp<-character()
          
          for(i in 1:wan_len){
            chr_pos<-which(genraw[,2]==wan[i,1])
            new_matrix<-genraw[chr_pos,]
            posi_pos<-which(new_matrix[,3]==wan[i,2])
            mark<-matrix(new_matrix[posi_pos,1],1,)
            marker<-rbind(marker,mark)
            sn<-matrix(new_matrix[posi_pos,4],1,)
            snp<-rbind(snp,sn)
          }
        }
        if(inputform==2){
          
          genRaw<-as.data.frame(genRaw)
          genraw<-genRaw[-1,1:4]
          
          wan_len<-dim(wan)[1]
          marker<-character()
          snp<-character()
          for(i in 1:wan_len){
            chr_pos<-which(genraw[,2]==wan[i,1])
            new_matrix<-genraw[chr_pos,]
            posi_pos<-which(new_matrix[,3]==wan[i,2])
            mark<-matrix(new_matrix[posi_pos,1],1,)
            marker<-rbind(marker,mark)
            sn<-matrix(new_matrix[posi_pos,4],1,)
            snp<-rbind(snp,sn)
          }
          
        }
        if(inputform==3){
          genRaw<-as.data.frame(genRaw)
          genraw<-genRaw[-1,c(1,3,4,12)]
          
          wan_len<-dim(wan)[1]
          marker<-character()
          snp<-character()
          for(i in 1:wan_len){
            chr_pos<-which(genraw[,2]==wan[i,1])
            new_matrix<-genraw[chr_pos,]
            posi_pos<-which(new_matrix[,3]==wan[i,2])
            mark<-matrix(new_matrix[posi_pos,1],1,)
            marker<-rbind(marker,mark)
            sn<-matrix(new_matrix[posi_pos,4],1,)
            snp<-rbind(snp,sn)
          }
        }
        
        wan<-cbind(marker,wan,maf,snp,vees,pees)
        
        tempwan <- wan
        lodscore1 <- as.numeric(tempwan[,5])
        log10P <- as.matrix(round(-log10(pchisq(lodscore1*4.605,1,lower.tail = F)),4))
        if(nrow(tempwan)>1){
          tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10]) 
        }else{
          tempwan1 <- cbind(t(as.matrix(tempwan[,1:5])),log10P,t(as.matrix(tempwan[,6:10])))  
        }
        wan <- tempwan1
        colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen (total)")
        wan<-as.data.frame(wan)
      }
    }
  }  
} 
parmsShow<-parmsShow[,-c(4,5,6,8,9)]
output<-list(result1=parmsShow,result2=wan)
return(output) 
}
}

