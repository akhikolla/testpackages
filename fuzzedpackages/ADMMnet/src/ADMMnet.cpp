// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X){
  int i, p=X.cols(), N=X.rows();
  Eigen::VectorXd mX(p), sdX(p);
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  return List::create(Named("x")=X, Named("sd")=sdX, Named("m")=mX);
}

/*****  Omega  *****/
// [[Rcpp::export]]
List OmegaC(Eigen::MatrixXd & Omega, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);
  
  //Omega.diagonal().setZero();
  Eigen::SparseMatrix<double> OmegaS=Omega.sparseView();
  
  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }
  
  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }
  
  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}

/*****  Sparse Omega  *****/
// [[Rcpp::export]]
List OmegaSC(Eigen::SparseMatrix<double> & OmegaS, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);
  
  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }
  
  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }
  
  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}





/////////////////////////////////
/////   Linear regression   /////
/////////////////////////////////

/*****  LM: Lambda path (max) inner product <xj,y> *****/
// [[Rcpp::export]]
double maxLambdaLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd wbeta, int N0){
  Eigen::VectorXd Li(N0);
  
  Li=(y.transpose()*X).cwiseAbs()/N0;   // <xj,y>/N0
  Li=Li.array()/wbeta.array()/alpha ;
  return Li.maxCoeff();
}


/*****  Used for CV trimming  *****/   //  
// [[Rcpp::export]]
Eigen::VectorXd cvTrimLmC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco, 
Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF) {
  int i, j;
  Eigen::VectorXd RSS, xbF=Eigen::VectorXd::Zero(NF); // xb=Eigen::VectorXd::Zero(N),
  
  if(nn2>0){
    RSS.setZero(nn2); //nn2= # of part of data
    
    if(nn==0){
      RSS(0)=yF.squaredNorm();
    }else{
      for(i=0;i<nn;i++){
        j=loco(i); //   index of nonzero beta
        xbF+=XF.col(j)*beta(i); // 
        RSS(i)=(yF-xbF).squaredNorm();
      }
    }
    if(nn2>nn){for(i=nn;i<nn2;i++){RSS(i)=RSS(nn-1);}}
  }else{
    RSS.setZero(1);
    RSS(0)=yF.squaredNorm();
  }
  
  return(RSS);
}


/*****  LM: Enet (L1+L2)  *****/
// [[Rcpp::export]]
List EnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
int isd, int p, int N0, double thresh, int maxit){ 
  
  int  i, j, it=0, il, iadd, ia=0; 
  double lambda2, zi, obj0, obj1, b0, db0, objQi, objQj, rss0, rss1, RSS0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd RSS(nlambda), RSQ(nlambda);
  double xr, dbMax, thresh2=1e-5;
  Eigen::VectorXd mX(p), sdX(p), di(p);
  
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();
  
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0);
    }
    lambda=lambda.array()*lambdaMax;
  }
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){  
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    local:
      while(1){
        ++it;
        
        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j); 
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update
        
        obj0=rss0/N0/2.0+objQj+objQi*lambda2/2.0;
        
        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}
      }//while
      
      iadd=0;
      for(i=0;i<p;++i){
        if(iactive(i)==0){
          di(i)=std::abs(y.dot(X.col(i))/N0);
          if(di(i)>lambda1(i)){
            active(ia)=i;iactive(i)=1;++ia;iadd=1;
          }
        }
      }
    if(iadd==1){goto local;}
    
    if (isd == 1) {
      Beta.col(il)=beta0; 
    } else {
      Beta.col(il)=beta0.array()/sdX.array();
    }
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
  return(List::create(Named("Beta")=Beta, Named("flag")=flag, Named("rsq")=RSQ,
  Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}


/*****  LM: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF){
  
  int i, j, it=0, il, iadd, ia=0; 
  double lambda2, zi, obj0, obj1, rss0, rss1, b0, db0, objQi, objQj, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda); // beta matrix for different lambdas
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXd RSSp(nlambda), RSS(nlambda), RSQ(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  double xr, dbMax, thresh2=1e-5;
  Eigen::VectorXd mX(p), sdX(p), di(p);
  
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha) 
    
      for(i=0;i<p;++i){
        if(iactive(i)==0){
          if(di(i)>lambda1(i)){
            active(ia)=i;iactive(i)=1;++ia;iadd=1;
          }
        }
      }
    
    //it=0;
    local:
      while(1){
        ++it;
        
        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j); 
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update
        
        obj0=rss0/N/2.0+objQj+objQi*lambda2/2.0;
        
        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}
      }//while
      
      iadd=0;
      for(i=0;i<p;++i){
        if(iactive(i)==0){
          di(i)=std::abs(y.dot(X.col(i))/N);
          if(di(i)>lambda1(i)){
            active(ia)=i;iactive(i)=1;++ia;iadd=1;
          }
        }
      }
    if(iadd==1){goto local;}
    
    Beta.col(il)=beta0.array()/sdX.array(); RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*(beta0(j)/sdX(j));}
    RSSp(il)=(yF-xbF).squaredNorm();
    
    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
  return(List::create(Named("Beta")=Beta, Named("flag")=flag,
  Named("RSS")=RSS, Named("rsq")=RSQ, Named("RSSp")=RSSp, Named("nlambda")=il));
}



/*****  LM: Network (L1+La)  *****/
// [[Rcpp::export]]
List NetLmC(Eigen::MatrixXd & X, Eigen::VectorXd & y, double alpha,
Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
int isd, int p, int N0, double thresh, int maxit){ 
  
  int i, j, ij, m, it=0, il, iadd, ia=0;
  double lambda2, zi, zi2, objQi=0.0, objQj, obj0, obj1, rss0, rss1, b0, db0, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N0);
  Eigen::VectorXd RSS(nlambda), RSQ(nlambda);
  double xr, dbMax, thresh2=1e-5, lambdaMax;
  Eigen::VectorXd mX(p), sdX(p), di(p);
  
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();
  
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0);
    }
    lambda=lambda.array()*lambdaMax;
  }
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    local:
      while(1){
        ++it;
        
        objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);} // Omega: w_kl/sqrt(d_k*d_l),L=I-Omega ; L=SLS (included sign of beta)
          }
          zi+=lambda2*zi2;
          
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1); 
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2); //  beta^T*L*beta
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2);
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2*zi2);
              beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update
        
        obj0=rss0/N0/2.0+objQj+objQi*lambda2/2.0;
        
        if(std::abs(rss1-rss0)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;} 
      }//while
      
      iadd=0;
      for(i=0;i<p;++i){
        if(iactive(i)==0){
           zi2=0.0;
           for(ij=0;ij<nadj(i);++ij){
              m=loc(ij, i);
              if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
           }
          di(i)=std::abs(y.dot(X.col(i))/N0+lambda2*zi2);
          if(di(i)>lambda1(i)){
            active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    if (isd == 1) {
      Beta.col(il)=beta0; 
    } else {
      Beta.col(il)=beta0.array()/sdX.array();
    }
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
  return(List::create(Named("Beta")=Beta, Named("flag")=flag, Named("rsq")=RSQ,
  Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  LM: Network (L1+La) cross-validation *****/
// [[Rcpp::export]]
List cvNetLmC(Eigen::MatrixXd & X, Eigen::VectorXd & y,double alpha,
Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF){ 
  
  int i, j, ij, m, it=0, il, iadd, ia=0;
  double lambda2, zi, zi2, objQi=0.0, objQj, obj0, obj1, rss0, rss1, b0, db0, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  Eigen::VectorXd RSS(nlambda), RSQ(nlambda), RSSp(nlambda);
  double xr, dbMax, thresh2=1e-5;
  Eigen::VectorXd mX(p), sdX(p), di(p);
  
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();
  
  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;
  
  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1-alpha);
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    local:
      while(1){
        ++it;
        
        objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);} // Omega: w_kl/sqrt(d_k*d_l),L=I-Omega ; L=SLS (included sign of beta)
          }
          zi+=lambda2*zi2;
          
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1); 
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2); //  beta^T*L*beta
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2);
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2*zi2);
              beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }
          
          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N/2.0+objQj+objQi*lambda2/2.0;
        
        if(std::abs(rss1-rss0)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;} 
      }//while
      
      iadd=0;
      for(i=0;i<p;++i){
        if(iactive(i)==0){
           zi2=0.0;
           for(ij=0;ij<nadj(i);++ij){
              m=loc(ij, i);
              if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
           }
          di(i)=std::abs(y.dot(X.col(i))/N+lambda2*zi2);
          if(di(i)>lambda1(i)){
            active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    Beta.col(il)=beta0.array()/sdX.array(); RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;
    
    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*(beta0(j)/sdX(j));}
    RSSp(il)=(yF-xbF).squaredNorm();
    
    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda
  
  exit:
  return(List::create(Named("Beta")=Beta, Named("flag")=flag, 
  Named("RSS")=RSS, Named("RSSp")=RSSp, Named("rsq")=RSQ, Named("nlambda")=il));
}




///////////////////
/////   Cox   /////
///////////////////



/*****  Cox: Lambda path (max)  *****/
// [[Rcpp::export]]
double maxLambdaCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, int N,
Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, double alpha, Eigen::VectorXd wbeta, int N0){
  int i, j, q;
  double denS=N, c1=0.0;
  Eigen::VectorXd Li(N), lli(N);
  
  for(i=0;i<n;i++){
    c1+=(nevent1(i)/denS);denS-=nevent(i);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){lli(j)=tevent(j)-c1;}
  }
  Li=(lli.transpose()*X).cwiseAbs()/N0;
  Li=Li.array()/wbeta.array()/alpha;
  return Li.maxCoeff();
}



/*****  Derivatives of log-pl of eta (1st&2nd order),  ties  *****/
void dletaCm(Eigen::VectorXd& exb, Eigen::VectorXd& tevent, int& N, 
Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1,
int& n, Eigen::VectorXd& pl1, Eigen::VectorXd& pl2, int& ifast, int& itwo){
  int i, j, q, ipl2=0;
  double denSi, c1=0.0, c2=0.0;
  Eigen::VectorXd denS(n);
  
  if(ifast==0 || itwo==1)goto two;
  denSi=exb.sum();
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denSi);c2+=(nevent1(i)/pow(denSi, 2));
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      denSi-=exb(j);
      pl1(j)=tevent(j)-exb(j)*c1;
      pl2(j)=exb(j)*(c1-exb(j)*c2);
      if(pl2(j)<=0.0)ipl2=1;
    }
  }
  if(ipl2==1){itwo=1;if(ifast==0){goto two;}}
  return;
  
  two:
  denSi=0.0;c1=0.0;c2=0.0;
  for(i=n-1;i>=0;--i){
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){denSi+=exb(j);}
    denS(i)=denSi;
  }
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denS(i));c2+=(nevent1(i)/pow(denS(i), 2));
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      pl1(j)=tevent(j)-exb(j)*c1;
      pl2(j)=exb(j)*(c1-exb(j)*c2);
    }
  }
}

/*****  Log-pl of eta,  ties  *****/
// [[Rcpp::export]]
double pletaCm(Eigen::VectorXd& xb, Eigen::VectorXd& exb, Eigen::VectorXi& nevent, 
Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n, int& ifast, int& itwo){
  int i, j, q, iSS=0;
  double ll=0.0, SSi;
  Eigen::VectorXd SS(n);

  if(ifast==0 || itwo==1)goto two;
  SSi=exb.sum();
  for(i=0;i<n;++i){
    if(SSi<=0.0)iSS=1;
    for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){ll+=xb(j);}
    ll-=nevent1(i)*log(SSi);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi-=exb(j);}
  }
  if(iSS==1){itwo=1;if(ifast==0){goto two;}}
  return(ll);
  
  two:
  ll=0.0;SSi=0.0;
  for(i=n-1;i>=0;--i){
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi+=exb(j);}
    SS(i)=SSi;
  }
  for(i=0;i<n;++i){
    for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){
      ll+=xb(j)-log(SS(i));
    }
  }
  return(ll);
}


/*****  Cox: Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvTrimCoxC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco, 
Eigen::MatrixXd XF, int NF, 
Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF, 
Eigen::MatrixXd X, int N, 
Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n, 
int ifast, int itwo){
  int i, j;
  double lli, lfi;
  Eigen::VectorXd cv, xb=Eigen::VectorXd::Zero(N), xbF=Eigen::VectorXd::Zero(NF);
  Eigen::VectorXd exb(N), exbF(NF);
  
  if(nn2>0){
    cv.setZero(nn2);
    
    if(nn==0){
      exb=(xb.array()).exp();
      lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      exbF=(xbF.array()).exp();
      lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
      cv(0)=lfi-lli;
    }else{
      for(i=0;i<nn;i++){
        j=loco(i);
        xb+=X.col(j)*beta(i);exb=(xb.array()).exp();
        lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        xbF+=XF.col(j)*beta(i);exbF=(xbF.array()).exp();
        lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
        cv(i)=lfi-lli;
      }
    }
    if(nn2>nn){for(i=nn;i<nn2;i++){cv(i)=cv(nn-1);}}
  }else{
    cv.setZero(1);
    
    exb=(xb.array()).exp();
    lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
    exbF=(xbF.array()).exp();
    lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
    cv(0)=lfi-lli;
  }
  
  return(cv);
}


/*****  Cox: Enet (L1+L2)  *****/
// [[Rcpp::export]]
List EnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, 
double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast){
  
  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
  double xi;
  
  for (i=0;i<p;++i) {
    xi=X.col(i).mean();
    X.col(i)=X.col(i).array()-xi;
  }
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQi=0.0;objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("nlambda")=il));
}


/*****  Cox: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, 
double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, 
int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){
  
  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
  double xi;
  
  for (i=0;i<p;++i) {
    xi=X.col(i).mean();
    X.col(i)=X.col(i).array()-xi;
  }
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQi=0.0;objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
          objQi+=pow(b0, 2);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
    
    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*beta0(j);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}



/*****  Cox: Network (L1+La)  *****/
// [[Rcpp::export]]
List NetCoxC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha, 
Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast){
  
  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        zi2=0.0;
        for(ij=0;ij<nadj(j);++ij){
          m=loc(ij, j)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
        }
        zi+=lambda2i*zi2;
        
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("nlambda")=il));
}


/*****  Cox: Network (L1+La)  cross-validation  *****/
// [[Rcpp::export]]
List cvNetCoxC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha, 
Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta, 
Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj, 
int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, 
int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF, 
int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){
  
  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
  
  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;
  
  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta;lambda2=lambda(il)*(1-alpha);
    lambda1i=N0*lambda1;lambda2i=N0*lambda2;
    
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }
    
    it=0;
    local:
    while(1){
      ++it;
      
      objQj=0.0;
      for(i=0;i<ia;++i){
        j=active(i);
        PLi2=pl2.dot(X.col(j).cwiseAbs2());
        zi=beta0(j)*PLi2+pl1.dot(X.col(j));
        zi2=0.0;
        for(ij=0;ij<nadj(j);++ij){
          m=loc(ij, j)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
        }
        zi+=lambda2i*zi2;
        
        if(zi>lambda1i(j)){
          b0=(zi-lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else if(zi<-lambda1i(j)){
          b0=(zi+lambda1i(j))/(lambda2i+PLi2);
          db0=beta0(j)-b0;
          objQi-=db0*(beta0(j)+b0-2*zi2);
          beta0(j)=b0;
          pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
          xb-=db0*X.col(j);
          objQj+=std::abs(b0)*lambda1(j);
        }else{
          b0=0.0;
          if(beta0(j)!=b0){
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
          }
        }
      }//for update
      
      ll1=ll0;obj1=obj0;
      exb=(xb.array()).exp();
      ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
      obj0=-ll0/N0+objQj+objQi*lambda2/2.0;
      
      if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
      if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
      if(obj0!=obj0){flag(il)=2;goto exit;}
      if(it>=maxit){flag(il)=1;goto exit;}
      
      dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
      if(ifast==1 && itwo==1)goto exit;
    }//while
    
    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}
    
    locbeta(il)=ll0;Betas.col(il)=beta0;
    
    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*beta0(j);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda
  
  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("flag")=flag, 
  Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}




