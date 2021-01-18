ISIS<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO){
pvalue<-svpal
svlod<-svmlod
inputform<-Genformat

if(is.null(psmatrix)){
  flagps<-1
}else{ 
  flagps<-0
}

y<-phe
ps<-psmatrix


if(is.null(svpal)==TRUE||is.null(svmlod)==TRUE){
  warning("Please set parameters!")
}

if((svpal<0)||(svpal>1))
{
  warning("Please input critical P-value between 0 and 1!")
}
if(svmlod<0)
{
  warning("Please input critical LOD score: > 0 !")
}
if(is.null(gen)==TRUE)
{
  warning("Please input correct genotypic dataset !","Warning",icon="warning")
  
}
if(is.null(y)==TRUE)
{
  warning("Please input correct phenotypic dataset !","Warning",icon="warning")
  
}

if((is.null(gen)==FALSE)&&(is.null(y)==FALSE)&&(ncol(gen)!=(nrow(y)+2)))
{
  warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset!","Error",icon="error")
  
}

if((is.null(gen)==FALSE)&&(is.null(y)==FALSE)&&((ncol(gen)==(nrow(y)+2)))&&(svpal>=0)&&(svpal<=1)&&(svmlod>=0))
{
wan<-NULL
result<-NULL

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

X<-t(gen)
set.seed(1)
X1<-X[3:nrow(X),]

sig<-seq(1:ncol(X1))
x<-data.frame(X1)
y<-as.matrix(y)
le<-length(sig)
xnew<-x[sig]

rm(x)
gc()

pval<-pvalue
y1<-matrix(nrow=nrow(y),ncol=ncol(y))

if (is.null(ps)==FALSE)
{
  ps1<-cbind(matrix(1,nrow=nrow(y)),ps)
  vhat<-solve(crossprod(ps1,ps1))%*%crossprod(ps1,y)
  vhat1<-vhat[-1]
  y1<-y-ps%*%vhat1
}else{
  y1<-y 
}
y1<-as.matrix(y1)

xxxq<-as.matrix(xnew)
rm(xnew)
gc()

mat<-vector()
matcor<-vector()

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

mm=foreach(i=1:le, .multicombine=TRUE, .combine = 'rbind')%dopar%
{
  if (var(xxxq[,i])>0){
    mat[i]<-cor.test(xxxq[,i],y1)$p.value
    matcor[i]<-abs(cor(xxxq[,i],y1))
    m<-c(mat[i],matcor[i])
  }else{
    mat[i]<-1
    matcor[i]<-0
    m<-c(mat[i],matcor[i])
  }
}
rownames(mm)<-NULL
mat<-mm[,1];matcor<-mm[,2]
stopCluster(cl)

rm(xxxq)
gc()


if(length(which(mat<pval))<=nrow(y)){
  ee<-as.vector(which(mat<pval))
}else{
  n1<-nrow(y1)-1
  ee<-as.vector(which(rank(matcor)>=(le-n1), arr.ind=T))
}

xxxnew<-X1[,sig[ee]]
yyy<-y1

set.seed(1)
cvfit1 <- ncvreg::cv.ncvreg(scale(xxxnew), yyy, family="gaussian",penalty="SCAD",gamma=3.7,warn=FALSE)


fit1 <- cvfit1$fit
obj11 <- (as.vector(fit1$beta[,cvfit1$min]))
obj1<-obj11[-1]
if (length(which(abs(obj1)!=0))!=0)
{
  sig1a<-which(abs(obj1)!=0)
  sig1b<-sig[-(sig[ee][sig1a])]
  yyy1<-y1-(X1[,sig[ee][sig1a]]%*%as.matrix(obj1[sig1a]))
  xxx1<-X1[,sig1b]
  mat1<-vector()
  for (i in 1:length(sig1b))
  {
    if (var(xxx1[,i])>0){
      
      mat1[i]<-abs(cor(xxx1[,i],yyy1))
    }else{
      mat1[i]<-0
      
    }
  }
  
  n2<-nrow(yyy1)-1
  ee1<-as.vector(which(rank(mat1)>=(ncol(xxx1)-n2), arr.ind=T))
  
  xxxnew1<-X1[,sig1b[ee1]]
  cvfit2 <- ncvreg::cv.ncvreg(scale(xxxnew1), yyy1, family="gaussian",penalty="SCAD",gamma=3.7,warn=FALSE)
  
  
  fit2 <- cvfit2$fit
  obj22 <- (as.vector(fit2$beta[,cvfit2$min]))
  
  rm(cvfit2,fit2,xxx1,xxxnew1)
  gc()
  
  obj2<-obj22[-1]
  sig1c<-sig1b[ee1][which(abs(obj2)!=0)]
  sigg<-sort(c(sig[ee][sig1a],sig1c))
  
}else{
  
  sigg<-sig[ee]
}
le1<-length(sigg)

if(le1!=0){
  
  ###########if result just have one column##############modified 2017.3.22###############
  if(le1==1){
    xxxnew11<-matrix(X1[,sigg],,1)
  }else{
    xxxnew11<-X1[,sigg]
  }
  
  z<-matrix()
  
  if (is.null(ps)==TRUE)
  {
    z<-matrix(1,nrow(X1),1)
  }else{
    z<-cbind(matrix(1,nrow(X1),1),ps) 
  }
  
  u1<-ebayes_EM(z,xxxnew11,y)
  obj3<-u1$u 
  result1<-matrix(0,ncol(X1)*1,ncol=1,nrow=ncol(X1))
  for (i in 1: le1)
  {
    result1[(sigg)[i],1]=obj3[i]
  }
  Res<- t(as.matrix((rowSums(result1)/ncol(result1))))
  Res1<-as.vector(Res)	
  le2<-length(which(abs(Res1)>1e-5))
  
  if(le2!=0){
    
    sig1<-which(abs(Res1)>1e-5)
    bbo<-matrix(0,le2,1)
    for (i in 1:le2){
      bbo[i,]=Res1[sig1[i]]
    }
    
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
    
    if(length(sig1)!= 0){
      
      if(length(sig1)==1){
        xxxx<-as.matrix(X1[,sig1])
        
      }else{
        xxxx<-X1[,sig1]
      }
      
      yn<-as.matrix(y)
      xxn<-z
      
      lod<-likelihood(xxn,xxxx,yn,bbo)
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
          qqq[j,2]=X[1,ii[j]]
          qqq[j,3]=X[2,ii[j]]
          qqq[j,4]=result1[ii[j],]
          
          qqq[j,5]=lod[j]
          qqq[j,6]=her[j]
        }
        id<-which(qqq[,5]==0)
        
        if(length(id)!=dim(qqq)[1]){
          
          if(length(id)!=0){
            qqq1<-qqq[-id,]
          }else{
            qqq1<-qqq
          }
          #######revised 2017 3.4##############################
          if(length(sig1)==1){
            xxmaf<-t(xxxx)
            xxmaf<-matrix(xxmaf,1,)
            result<-matrix(qqq1[,-1],1,)
          }else{
            xxmaf<-t(xxxx)
            result<-as.matrix(qqq1[,-1])
          }
          
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
            
            vees<-matrix("",nrow = nrow(result),1)
            pees<-matrix("",nrow = nrow(result),1)
            pees[1,1]<-pee
            vees[1,1]<-vee
            result<-as.matrix(qqq1[,-1])
            result<-result
            temp<-as.matrix(result[,3:5])
            temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
            temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
            wan<-cbind(result[,1:2],temp)
            
          }else{
            pees<-as.matrix(pee)
            vees<-as.matrix(vee) 
            result<-t(as.matrix(qqq1[,-1]))
            result<-result
            temp<-t(as.matrix(result[,3:5]))
            temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
            temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
            wan<-cbind(t(as.matrix(result[,1:2])),temp)
            
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
        }
      }
    }
  }  
}
output<-list(result=wan)
}
return(output)
#}
}