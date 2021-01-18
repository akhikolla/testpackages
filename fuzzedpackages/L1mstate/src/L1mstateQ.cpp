// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;
//typedef Eigen::SparseVector<double> SpVec;
//typedef Eigen::SparseVector<double> SpVecf;
//typedef Eigen::SparseVector<int> SpVeci;
//typedef SpVec::InnerIterator InIterVec;
//typedef Eigen::SparseMatrix<double> SpMat;


/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleQ(Eigen::MatrixXd X){
  Eigen::VectorXd mX=X.colwise().mean(), sdX(X.cols()), sdXi;
  X.rowwise()-=mX.transpose();
  sdX=X.colwise().norm()/sqrt((double)X.rows());
  sdXi=1.0/sdX.array();
  X=X*sdXi.asDiagonal();
  return List::create(Named("x")=X, Named("sd")=sdX);
}


/*****  Log-pl of eta,  ties  *****/
// [[Rcpp::export]]
double pletaQ(Eigen::VectorXd& xb, Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n){
  int i, j, q;
  double ll=0.0, SSi=0.0;
  
  if(loc1.size()==0){
    return 0.0;
  }else{
    Eigen::VectorXd exb = (xb.array()).exp();
    SSi=exb.sum();
    for(i=0;i<n;++i){
      for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){ll+=xb(j);}
      ll-=nevent1(i)*log(SSi);
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi-=exb(j);}
    }
    return(ll);
  }
}


/*****  1st order derivatives of log-pl of eta,  ties *****/
/*****  1st order derivatives of negavitve-log-pl of eta = -pl1*****/
// [[Rcpp::export]]
Eigen::VectorXd d1Q(Eigen::VectorXd& xb, Eigen::VectorXd& tevent, int& N, Eigen::VectorXi& nevent, 
                    Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n){
  int i, j, q;
  double denSi, c1=0.0;
  Eigen::VectorXd pl1(N);
  Eigen::VectorXd exb = (xb.array()).exp();
  denSi=exb.sum();
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denSi);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      denSi-=exb(j);
      pl1(j)=tevent(j)-exb(j)*c1;
    }
  }
  return(pl1);
}


/*****  2nd order derivatives of log-pl of eta,  ties *****/
/*****  2nd order derivatives of negavitve-log-pl of eta = -pl2*****/
// [[Rcpp::export]]
Eigen::VectorXd d2Q(Eigen::VectorXd& xb, Eigen::VectorXd& tevent, int& N, Eigen::VectorXi& nevent, 
                    Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n){
  int i, j, q;
  double denSi, c1=0.0, c2=0.0;
  Eigen::VectorXd pl2(N);
  Eigen::VectorXd exb = (xb.array()).exp();
  denSi=exb.sum();
  for(i=0;i<n;++i){
    c1+=(nevent1(i)/denSi);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
      denSi-=exb(j);
      pl2(j)=exb(j)*(c1-exb(j)*c2);
    }
  }
  return(pl2);
}


/***** Lambda path (max)  *****/
// [[Rcpp::export]]
double max_lambdaQ(List& x){
  int Q = x.size();
  int q, ia=0, N, n, k;
  Eigen::VectorXd L(Q), tevent, s(Q), Li;
  List xq;
  Eigen::VectorXi nevent, nevent1, loc1;
  Eigen::MatrixXd xq0;
  
  for(q=0; q<Q; ++q){
    xq = x[q]; loc1 = xq[6];
    if(loc1.size()==0) {L(q)=0.0;}
    else {s(ia)=q;++ia;}
  }
  for(k=0; k<ia; ++k){
    q = s(k);
    xq = x[q]; xq0 = xq[0]; tevent = xq[2]; N = xq[3]; nevent = xq[4]; nevent1 = xq[5];loc1 = xq[6]; n = xq[7];
    Eigen::VectorXd xb = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd pl1 = Eigen::VectorXd::Zero(N);
    pl1 = d1Q(xb, tevent, N, nevent, nevent1, loc1, n);
    
    Li=(pl1.transpose()*xq0).cwiseAbs();
    L(q)=Li.maxCoeff()/N;
  }
  return (L.maxCoeff());
}


/*****  LASSO  *****/
// [[Rcpp::export]]
List l1msQ(List& x, Eigen::VectorXd lambda, int nlambda, int p, double thresh, int maxit){
  //number of transitions
  int Q = x.size();
  //information for each transition
  List xq(8);
  Eigen::MatrixXd xq0;
  Eigen::MatrixXd X;
  Eigen::VectorXd tevent;
  Eigen::VectorXi nevent, nevent1, loc1;
  int Nq, N, n;
  //intermediate quantities
  //beta from previous iteration and current beta
  Eigen::MatrixXd beta0Q = Eigen::MatrixXd::Zero(Q,p);
  Eigen::MatrixXd betacQ = Eigen::MatrixXd::Zero(Q,p);
  Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd betac = Eigen::VectorXd::Zero(p);
  double lambda1;
  Eigen::VectorXd lambda1i = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd W = Eigen::VectorXd::Constant(p,1.0);
  //log-likelihood values
  Eigen::VectorXd llQ = Eigen::VectorXd::Zero(Q);
  //eta and first-and second-order derivatives
  Eigen::VectorXd xb, exb, pl1, pl2;
  //active set
  Eigen::VectorXi active;
  Eigen::VectorXi iactive;
  int iadd;
  //others
  int il, q, it, j;
  double ll0=0.0, obj0=0.0, objQj, b0, db0, ll1, obj1, PLi, PLi2, zi;
  
  //output
  Eigen::VectorXd locbeta = Eigen::VectorXd::Zero(nlambda);
  List Betas(nlambda);
  
  //Initialization
  // locbeta(0)
  for(q=0;q<Q;++q){
    xq = x[q]; 
    X = xq[0]; tevent = xq[2]; N = xq[3]; nevent = xq[4]; nevent1 = xq[5]; loc1 = xq[6]; n = xq[7];
    if(loc1.size()==0){
      llQ(q) = 0.0;
    }else{
      xb = Eigen::VectorXd::Zero(N);
      llQ(q) = pletaQ(xb,nevent,nevent1,loc1,n);
    }
  }
  locbeta(0) = llQ.sum();
  Betas[0] = Eigen::MatrixXd::Zero(Q,p);
  
  //lambda path
  for(il = 1; il<nlambda; ++il){
    //lambda value
    lambda1 = lambda(il);
    //assign beta0Q 
    beta0Q = Betas[il-1];
    betacQ = Eigen::MatrixXd::Zero(Q,p);
    //q-loop
    q=0;
    label:
      while(q>=0 && q<Q){
        //information for q-transition 
        xq = x[q]; X = xq[0]; Nq=xq[1]; tevent = xq[2]; N = xq[3]; nevent = xq[4]; nevent1 = xq[5];loc1 = xq[6]; n = xq[7];
        
        //check whether or not q-transition happens
        if(loc1.size()==0){
          //results of q-transition
          for(int i=0; i<p; ++i){
            betacQ(q,i) = 0.0;
          }
          llQ(q) = 0.0;
          
          goto label1;
        }else{
          lambda1i = Eigen::VectorXd::Zero(p);
          lambda1i = Nq*lambda1*W;
          //assign beta0
          betac = Eigen::VectorXd::Zero(p);
          for(int i=0; i<p; i++){beta0(i) = beta0Q(q,i);}
          //active set
          active=Eigen::VectorXi::Zero(p);
          for(int i=0;i<p;++i){
            if(active(i)==0){
              if(beta0(i) != 0){active(i) = 1;}
            }
          }
          
          
          it=0;
          //while loop
          local:
            while(1){
              it++;
              
              xb = Eigen::VectorXd::Zero(N);
              for(int i=0; i<p; i++){xb += (X.col(i))*beta0(i);}
              //calculate first- and second-order derivatives pl1, pl2
              pl1 = d1Q(xb, tevent, N, nevent, nevent1, loc1, n);
              pl2 = d2Q(xb, tevent, N, nevent, nevent1, loc1, n);
              
              objQj=0.0;
              for(j=0;j<p;++j){
                if(active(j)){
                  PLi2=pl2.dot(X.col(j).cwiseAbs2());
                  zi=beta0(j)*PLi2+pl1.dot(X.col(j));
                  if(zi>lambda1i(j)){
                    b0=(zi-lambda1i(j))/(PLi2);
                    db0=beta0(j)-b0;betac(j)=b0;
                    pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
                    xb-=db0*X.col(j);
                    objQj+=std::abs(b0)*lambda1;
                  }else if(zi<-lambda1i(j)){
                    b0=(zi+lambda1i(j))/(PLi2);
                    db0=beta0(j)-b0;betac(j)=b0;
                    pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
                    xb-=db0*X.col(j);
                    objQj+=std::abs(b0)*lambda1;
                  }else{
                    b0=0.0;
                    if(beta0(j)!=b0){
                      db0=beta0(j)-b0;betac(j)=b0;
                      pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
                      xb-=db0*X.col(j);
                    }
                  }
                }
              }//for update
              
              ll1=ll0;obj1=obj0;
              exb=(xb.array()).exp();
              ll0=pletaQ(xb, nevent, nevent1, loc1, n);
              obj0=-ll0/Nq+objQj;
              //check convergence
              for(int i=0; i<p; i++) beta0(i)=betac(i);
              if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){break;}
              if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){break;}
              if(obj0!=obj0){goto label2;}
              if(it>=maxit){goto label2;}
              
            }//end-while-loop
            
            //scan for active  
            iadd=0;
          for(int i=0;i<p;++i){
            if(active(i)==0){
              PLi=pl1.dot(X.col(i));
              if(std::abs(PLi)>lambda1i(i)){
                active(i)=1;iadd=1;
              }
            }
          }
          if(iadd==0){
            goto label2;
          }else{
            for(int j=0; j<p; j++){beta0(j) = betac(j);}
            goto local;
          }
          
          label2:
            //results of q-transition
            for(int i=0; i<p; i++){betacQ(q,i) = betac(i);}
            xb = Eigen::VectorXd::Zero(N);
          for(int i=0; i<p; i++){xb += (X.col(i))*betac(i);}
          llQ(q) = pletaQ(xb,nevent,nevent1,loc1,n); 
        }
        
        label1:
          if(q==(Q-1)){
            goto exit;
          }else{
            q=q+1; goto label;
          }
      }//end-q(transition)-loop
      
      exit:
        //results of all Q transitions  
        locbeta(il)=llQ.sum(); Betas[il]=betacQ;
  }//end-lambda-path
  
  return(List::create(Named("Beta")=Betas, Named("ll")=locbeta, Named("nlambda")=il));
}

/*****  LASSO  *****/
/*****  cross-validation PL  *****/
// [[Rcpp::export]]
List cvl1msQ(List& x, Eigen::VectorXd lambda, int nlambda, int p, double thresh, int maxit, List& xF){
  //number of transitions
  int Q = x.size();
  //information for each transition
  List xq(8);
  Eigen::MatrixXd xq0;
  Eigen::MatrixXd X;
  Eigen::VectorXd tevent;
  Eigen::VectorXi nevent, nevent1, loc1;
  int Nq, N, n;
  //information for each transition of xF
  List xqF(12);
  Eigen::MatrixXd XF;
  Eigen::VectorXd teventF;
  Eigen::VectorXi neventF, nevent1F, loc1F;
  int NF, nF;
  //intermediate quantities
  //beta from previous iteration and current beta
  Eigen::MatrixXd beta0Q = Eigen::MatrixXd::Zero(Q,p);
  Eigen::MatrixXd betacQ = Eigen::MatrixXd::Zero(Q,p);
  Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd betac = Eigen::VectorXd::Zero(p);
  double lambda1;
  Eigen::VectorXd lambda1i = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd W = Eigen::VectorXd::Constant(p,1.0);
  //log-likelihood values
  Eigen::VectorXd llQ = Eigen::VectorXd::Zero(Q);
  Eigen::VectorXd llF = Eigen::VectorXd::Zero(Q);
  //eta and first-and second-order derivatives
  Eigen::VectorXd xb, exb, pl1, pl2;
  Eigen::VectorXd xbqF = Eigen::VectorXd::Zero(NF);
  //active set
  Eigen::VectorXi active;
  Eigen::VectorXi iactive;
  int iadd;
  //others
  int il, q, it, j;
  double ll0=0.0, obj0=0.0, objQj, b0, db0, ll1, obj1, PLi, PLi2, zi;
  
  //output
  Eigen::VectorXd locbeta = Eigen::VectorXd::Zero(nlambda);
  Eigen::VectorXd locbetaF = Eigen::VectorXd::Zero(nlambda);
  List Betas(nlambda);
  
  //Initialization
  // locbeta(0)
  for(q=0;q<Q;++q){
    xq = x[q]; 
    X = xq[0]; tevent = xq[2]; N = xq[3]; nevent = xq[4]; nevent1 = xq[5]; loc1 = xq[6]; n = xq[7];
    if(loc1.size()==0){
      llQ(q) = 0.0;
    }else{
      xb = Eigen::VectorXd::Zero(N);
      llQ(q) = pletaQ(xb,nevent,nevent1,loc1,n);
    }
  }
  locbeta(0) = llQ.sum();
  for(int k=0; k<Q; ++k){
    xqF = xF[k]; XF = xqF[0];
    NF = xqF[3]; neventF = xqF[4]; nevent1F = xqF[5];loc1F = xqF[6]; nF = xqF[7];
    if(loc1F.size()==0){
      llF(q) = 0.0;
    }else{
      xbqF = Eigen::VectorXd::Zero(NF);
      llF(k) = pletaQ(xbqF, neventF, nevent1F, loc1F, nF);
    }
  }
  locbetaF(0) = llF.sum();
  Betas[0] = Eigen::MatrixXd::Zero(Q,p);
  
  //lambda path
  for(il = 1; il<nlambda; ++il){
    //lambda value
    lambda1 = lambda(il);
    //assign beta0Q 
    beta0Q = Betas[il-1];
    betacQ = Eigen::MatrixXd::Zero(Q,p);
    //q-loop
    q=0;
    label:
      while(q>=0 && q<Q){
        //information for q-transition 
        xq = x[q]; X = xq[0]; Nq=xq[1]; tevent = xq[2]; N = xq[3]; nevent = xq[4]; nevent1 = xq[5];loc1 = xq[6]; n = xq[7];
        
        //check whether or not q-transition happens
        if(loc1.size()==0){
          //results of q-transition
          for(int i=0; i<p; ++i){
            betacQ(q,i) = 0.0;
          }
          llQ(q) = 0.0;
          
          goto label1;
        }else{
          lambda1i = Eigen::VectorXd::Zero(p);
          lambda1i = Nq*lambda1*W;
          //assign beta0
          betac = Eigen::VectorXd::Zero(p);
          for(int i=0; i<p; i++){beta0(i) = beta0Q(q,i);}
          //active set
          active=Eigen::VectorXi::Zero(p);
          for(int i=0;i<p;++i){
            if(active(i)==0){
              if(beta0(i) != 0){active(i) = 1;}
            }
          }
          
          it=0;
          //while loop
          local:
            while(1){
              it++;
              
              xb = Eigen::VectorXd::Zero(N);
              for(int i=0; i<p; i++){xb += (X.col(i))*beta0(i);}
              //calculate first- and second-order derivatives pl1, pl2
              pl1 = d1Q(xb, tevent, N, nevent, nevent1, loc1, n);
              pl2 = d2Q(xb, tevent, N, nevent, nevent1, loc1, n);
              
              objQj=0.0;
              for(j=0;j<p;++j){
                if(active(j)){
                  PLi2=pl2.dot(X.col(j).cwiseAbs2());
                  zi=beta0(j)*PLi2+pl1.dot(X.col(j));
                  if(zi>lambda1i(j)){
                    b0=(zi-lambda1i(j))/(PLi2);
                    db0=beta0(j)-b0;betac(j)=b0;
                    pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
                    xb-=db0*X.col(j);
                    objQj+=std::abs(b0)*lambda1;
                  }else if(zi<-lambda1i(j)){
                    b0=(zi+lambda1i(j))/(PLi2);
                    db0=beta0(j)-b0;betac(j)=b0;
                    pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
                    xb-=db0*X.col(j);
                    objQj+=std::abs(b0)*lambda1;
                  }else{
                    b0=0.0;
                    if(beta0(j)!=b0){
                      db0=beta0(j)-b0;betac(j)=b0;
                      pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
                      xb-=db0*X.col(j);
                    }
                  }
                }
              }//for update
              
              ll1=ll0;obj1=obj0;
              exb=(xb.array()).exp();
              ll0=pletaQ(xb, nevent, nevent1, loc1, n);
              obj0=-ll0/Nq+objQj;
              //check convergence
              for(int i=0; i<p; i++) beta0(i)=betac(i);
              if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){break;}
              if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){break;}
              if(obj0!=obj0){goto label2;}
              if(it>=maxit){goto label2;}
              
              //scan for active    
            }//end-while-loop
            
            iadd=0;
          for(int i=0;i<p;++i){
            if(active(i)==0){
              PLi=pl1.dot(X.col(i));
              if(std::abs(PLi)>lambda1i(i)){
                active(i)=1;iadd=1;
              }
            }
          }
          if(iadd==0){
            goto label2;
          }else{
            for(int j=0; j<p; j++){beta0(j) = betac(j);}
            goto local;
          }
          
          label2:
            //results of q-transition
            for(int i=0; i<p; i++){betacQ(q,i) = betac(i);}
            xb = Eigen::VectorXd::Zero(N);
          for(int i=0; i<p; i++){xb += (X.col(i))*betac(i);}
          llQ(q) = pletaQ(xb,nevent,nevent1,loc1,n); 
        }
        
        label1:
          if(q==(Q-1)){
            goto exit;
          }else{
            q=q+1; goto label;
          }
      }//end-q(transition)-loop
      
      exit:
        //results of all Q transitions  
        locbeta(il)=llQ.sum(); Betas[il]=betacQ;
        
        //compute the cross-validated log-likelihood of the whole dataset
        for(int k=0; k<Q; ++k){
          xqF = xF[k]; XF = xqF[0];
          NF = xqF[3]; neventF = xqF[4]; nevent1F = xqF[5];loc1F = xqF[6]; nF = xqF[7];
          
          xbqF = Eigen::VectorXd::Zero(NF);
          //Update eta of the whole dataset
          for(int i=0; i<p; i++){
            xbqF += XF.col(i)*betacQ(k,i);
          }
          
          //log-likelihood value
          llF(k) = pletaQ(xbqF, neventF, nevent1F, loc1F, nF);
        }
        locbetaF(il) = llF.sum();
        
  }//end-lambda-path
  
  return(List::create(Named("beta")=Betas, Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=nlambda));
}
