// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;


//////////////////////////////////////////////////////
/////   2-Paramter Logistic Model with L2-Norm   /////
//////////////////////////////////////////////////////

// [[Rcpp::export]]
List P2L2C(Eigen::MatrixXd X, Eigen::VectorXd y,
           Eigen::ArrayXd ldS, double ltp, Eigen::VectorXd lambda, Eigen::ArrayXd wldS, 
           double thresh, int maxit, double threshP, double threshB){
  
  int  i, j, i2, it=0, il, nld=ldS.size(), nlambda=lambda.size(), p=X.cols(), N0=X.rows();
  double zi, b0, db0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1=Eigen::VectorXd::Zero(p), lambda2(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  
  double dbMax;
  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);
  
  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;
  
  
  double ld0=0.0, ld02=0.0;
  double wld0=wldS.sum();
  
  
  //  Initial values
  mX(0)=0.0; sdX(0)=1.0;
  X2.col(0)=X2.col(0).array()+1.0;
  
  if (p>1) {
    i=1;
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
   
    for (i=0; i<nld; ++i) {
      ldS(i)=(ldS(i)-mX(1))/sdX(1);
    }
  }
  
  for (i=0; i<nld; ++i) {
    ld0+=ldS(i)*wldS(i);
    ld02+=ldS(i)*ldS(i)*wldS(i);
  }
  
  
  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  
  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  
  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();
  
  
  for (il=0; il<nlambda; ++il) {
    
    it=0;
    while(1){
      ++it;
      
      W=p0*(1.0-p0);
      Z=(y.array()-p0)/W.array();
      
      dbMax=0.0;
      for (j=0; j<p; ++j) {
        
        wi2=W.dot(X2.col(j))/N0;
        xr=(Z.array()*W.array()).matrix().dot(X.col(j));
        zi=xr/N0+wi2*beta0(j);
        
        if (j==0) {
          zi-=lambda(il)*(beta0(1)*ld0-ltp*wld0);
          b0=zi/(lambda(il)*wld0+wi2);
          
          if (b0 > threshB) {
            b0=threshB;
          } else if (b0 < (-threshB)) {
            b0=-threshB;
          }
          
        } else {
          zi-=lambda(il)*(beta0(0)-ltp)*ld0;
          b0=zi/(lambda(il)*ld02+wi2);
          
          if (b0 > threshB) {
            b0=threshB;
          } else if (b0 < 0.0) {
            b0=0.001;
          }
          
        }
        
        db0=b0-beta0(j); beta0(j)=b0;
        
        Z-=db0*X.col(j);
        xb+=db0*X.col(j).array();
        
        dbMax=std::max(dbMax, pow(db0, 2));
      }//for update
      
      
      for (i2=0; i2<N0; ++i2) {
        p0(i2)=1.0/(1.0+exp(-xb(i2)));
        if (p0(i2) < threshP) {
          p0(i2)=threshP;
        } else if (p0(i2) > (1.0-threshP)) {
          p0(i2)=1.0-threshP;
        }
      }
      
      if(dbMax<thresh){flag(il)=0; break;}
      if(it>=maxit){flag(il)=1; break;
      // goto exit;
      }
    }//while
    
    Z=(y.array()-p0);
    
    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();
    
    BetaSTD.col(il)=beta0; 
    Beta.col(il)=beta0.array()/sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));
    
  }//for lambda
  
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  Cross-Validation  *****/
// [[Rcpp::export]]
List cvP2L2C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::MatrixXd XF, Eigen::VectorXd yF,
             Eigen::ArrayXd ldS, double ltp, Eigen::VectorXd lambda, Eigen::ArrayXd wldS, 
             double thresh, int maxit, double threshP, double threshB){
  
  int  i, j, i2, it=0, il, nld=ldS.size(), nlambda=lambda.size(), p=X.cols(), N0=X.rows(), NF=XF.rows();
  double zi, b0, db0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;
  
  Eigen::VectorXd mX(p), sdX(p), di(p);
  
  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);
  
  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;
  
  Eigen::ArrayXd xbF(NF), pF0(NF), LLF(nlambda);
  double llF0;
  
  double ld0=0.0, ld02=0.0;
  double wld0=wldS.sum();
  
  
  //  Initial values
  mX(0)=0.0; sdX(0)=1.0;
  X2.col(0)=X2.col(0).array()+1.0;
  
  if (p>1) {
    i=1;
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
    
    for (i=0; i<nld; ++i) {
      ldS(i)=(ldS(i)-mX(1))/sdX(1);
    }
  }
  
  
  for (i=0; i<nld; ++i) {
    ld0+=ldS(i)*wldS(i);
    ld02+=ldS(i)*ldS(i)*wldS(i);
  }
  
  
  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  
  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  
  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();
  
  
  for(il=0;il<nlambda;++il){

    it=0;
    while(1){
      ++it;
      
      W=p0*(1.0-p0);
      Z=(y.array()-p0)/W.array();
      
      dbMax=0.0; 
      for(j=0;j<p;++j){
        
        wi2=W.dot(X2.col(j))/N0;
        xr=(Z.array()*W.array()).matrix().dot(X.col(j));
        zi=xr/N0+wi2*beta0(j);
        
        if (j==0) {
          zi-=lambda(il)*(beta0(1)*ld0-ltp*wld0);
          b0=zi/(lambda(il)*wld0+wi2);
          
          if (b0 > threshB) {
            b0=threshB;
          } else if (b0 < (-threshB)) {
            b0=-threshB;
          }
          
        } else {
          zi-=lambda(il)*(beta0(0)-ltp)*ld0;
          b0=zi/(lambda(il)*ld02+wi2);
          
          if (b0 > threshB) {
            b0=threshB;
          } else if (b0 < 0.0) {
            b0=0.001;
          }
          
        }
        
        db0=b0-beta0(j); beta0(j)=b0;
        
        Z-=db0*X.col(j);
        xb+=db0*X.col(j).array();
        
        dbMax=std::max(dbMax, pow(db0, 2));
      }//for update
      
      for (i2=0; i2<N0; ++i2) {
        p0(i2)=1.0/(1.0+exp(-xb(i2)));
        if (p0(i2) < threshP) {
          p0(i2)=threshP;
        } else if (p0(i2) > (1.0-threshP)) {
          p0(i2)=1.0-threshP;
        }
      }
      
      if(dbMax<thresh){flag(il)=0; break;}
      if(it>=maxit){flag(il)=1; break;
      // goto exit;
      }
    }//while
    
    Z=(y.array()-p0);
    
    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();
    
    BetaSTD.col(il)=beta0; 
    Beta.col(il)=beta0.array()/sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));
    
    
    // Predict Deviance
    xbF=XF*Beta.col(il);
    for (i2=0; i2<NF; ++i2) {
      pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
      if (pF0(i2) < threshP) {
        pF0(i2)=threshP;
      } else if (pF0(i2) > (1.0-threshP)) {
        pF0(i2)=1.0-threshP;
      }
    }
    LLF(il)=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();
    
  }//for lambda
  
  xbF=log(yF.mean()/(1.0-yF.mean()))*XF.col(0).array();
  for (i2=0; i2<NF; ++i2) {
    pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
    if (pF0(i2) < threshP) {
      pF0(i2)=threshP;
    } else if (pF0(i2) > (1.0-threshP)) {
      pF0(i2)=1.0-threshP;
    }
  }
  llF0=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();
  
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                            Named("LL")=LL.head(il), Named("ll0")=ll0,
                            Named("LLF")=LLF.head(il), Named("llF0")=llF0,
                            Named("lambda")=lambda, Named("nlambda")=il));
}




