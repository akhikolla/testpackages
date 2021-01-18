
#ifndef xxCOMMON_H
#define xxCOMMON_H

//#include <RcppArmadillo.h>

double CrossEntropyVal(arma::vec yVec, arma::vec tVec);
arma::mat CrossEntropy(arma::mat y, arma::mat t);
arma::mat NewCrossEntropy(arma::mat y, arma::mat t);
arma::mat Softmax(arma::mat X);
arma::mat Masking(arma::mat X, double cut_off);
arma::mat AlphaMasking(arma::mat X, double cut_off, double alpha);
arma::vec BatchNorm(arma::vec x);
arma::vec BackwardBatchNorm(arma::vec x, arma::vec dOut);

String Num2ActiveStr(int i);
void MakeStrVec(arma::vec nstrVec, String* strVec);
double GetAccuracy(arma::mat X, arma::mat t_Mat);
IntegerVector RandInts(int nMany, int ceiling);
int GetRemains(int a, int b);
arma::mat Con2OneHotEncoding(arma::mat X);


arma::mat fi(arma::mat v, String strDist);
arma::mat dfi(arma::mat v, arma::mat Out, arma::mat xdOut, String strDist);


arma::mat dfi(arma::mat v, arma::mat Out, arma::mat xdOut, String strDist){
  
  int q = xdOut.n_rows;           ///// V:q2xn  Out:qxn   dOut:q2xn  _dOut:qxn  
  int n = xdOut.n_cols;
  int q2 = q+2;
  
  arma::mat dOut(q2, n);
  
  double del=1e-5;
  
  double val1=0;
  double val2=0;
  double val3=0;
  
  double mu=0;
  double sig=0;                   ///// V:q2xn  Out:qxn   dOut:q2xn  _dOut:qxn  
  
  arma::vec x(q);
  
  for(int j=1;j<=n;j++){
    mu = v(0,j-1);
    sig = v(1,j-1);
    dOut(0,j-1) = accu(xdOut.col(j-1) );
    dOut(1,j-1) = accu( ( (Out.col(j-1) - mu)/(sig+del) )  % xdOut.col(j-1)  );
    
    if(strDist==strExponential){
      dOut(0,j-1) = accu( -log(1- v.submat(2,(j-1),(q2-1),(j-1)) ) % xdOut.col(j-1) );
    }else{
      mu = v(0,j-1);
      sig = v(1,j-1);
      dOut(0,j-1) = accu(xdOut.col(j-1) );
      dOut(1,j-1) = accu( ( (Out.col(j-1) - mu)/(sig+del) )  % xdOut.col(j-1)  );
      
    }
    
    for(int i=3;i<=q2;i++){
      val1 = Out(i-1-2,j-1);
      
      if(strDist == strNormal){
        val2 = R::dnorm(val1, mu, sig, 0); 
        val3 = 1/(val2+del);  
      }else if(strDist == strLogistic){
        val2 = R::dlogis(val1, mu, sig, 0);
        val3 = 1/(val2+del);  
      }else if(strDist == strCauchy){
        val2 = R::dcauchy(val1, mu, sig, 0);
        val3 = 1/(val2+del);  
      }else if(strDist == strExponential){
        val2 = R::dexp(val1, 1/mu, 0);
        val3 = 1/(val2+del);    
      }else{
        val2 = R::dnorm(val1, mu, sig, 0) ;
        val3 = 1/(val2+del);  
      }
      
      
      dOut(i-1,j-1) = val3 * xdOut(i-2-1,j-1); 
    }
    
  }
  
  return dOut;  
  
  
}


arma::mat fi(arma::mat v, String strDist){
  
  int q2 = v.n_rows;               /// V:q2xn 
  int q = q2-2;
  int n =  v.n_cols;
  double u=0;
  double mu=0;
  double sig=0;
  
  arma::mat Out(q,n);              /// Out :qxn   V-> FInv ->Out
  
  for(int j=1;j<=n;j++){
    mu = v(0,j-1);
    sig = v(1,j-1);
    
    for(int i=1;i<=q;i++){
      u = v(i+2-1, j-1);
      if(strDist == strNormal){
        Out(i-1, j-1) = R::qnorm(u, mu, sig, 1,0); 
      }else if(strDist == strLogistic){
        Out(i-1, j-1) = R::qlogis(u, mu, sig, 1,0); 
      }else if(strDist == strCauchy){
        Out(i-1, j-1) = R::qcauchy(u, mu, sig, 1,0); 
      }else if(strDist == strExponential){
        Out(i-1, j-1) = R::qexp(u, 1/mu, 1,0); 
      }else{
        Out(i-1, j-1) = R::qnorm(u, mu, sig, 1,0) ;
      }
      
    }
    
  }
  
  return Out;
  
  
  
}






arma::mat Con2OneHotEncoding(arma::mat X){
  
  int p=X.n_rows;
  int n=X.n_cols;
  
  arma::mat out(p,n);
  out.zeros();
  
  int nIndex=0;
  arma::vec x(p);
  for(int i=1;i<=n;i++){
    x = X.col(i-1);
    nIndex = x.index_max();
    out(nIndex, i-1) = 1;
  }
  
  return out;

}






IntegerVector RandInts(int nMany, int ceiling) {
  
  bool bMode = FALSE;
  if(nMany > ceiling){
    bMode = TRUE;
  }
  IntegerVector results(nMany) ;
  
  IntegerVector frame = seq_len(ceiling) ;
  
  IntegerVector candidate(nMany) ;
  int maxx=ceiling+1;
  
  while (maxx > ceiling) {
    
    candidate = RcppArmadillo::sample(frame, nMany, bMode, NumericVector::create() ) ;
    
    maxx = max(candidate);
    results = candidate;
    
  }
  
  return results;
}




int GetRemains(int a, int b){
  
  int quo = a/b;
  int remain = a-(b*quo);
  return remain;
}




double GetAccuracy(arma::mat X, arma::mat t_Mat){
  
  int p=X.n_rows;
  int n=X.n_cols;
  
  int nIndex=0;
  arma::vec x(p);
  double dAccuracy = 0;
  for(int i=1;i<=n;i++){
    x = X.col(i-1);
    nIndex = x.index_max();
    if( t_Mat(nIndex, i-1) == 1 ){
      dAccuracy += 1;
    }
  }
  
  double out = dAccuracy/n;
  return out;


}


String Num2ActiveStr(int i){
  
  String ans("");
  
  if(i==1){
    return strSigmoid;
  }else if(i==2){
    return strRelu;
  }else if(i==3){
    return strLeakyRelu;
  }else if(i==4){
    return strTanH;
  }else if(i==5){
    return strArcTan;
  }else if(i==6){
    return strArcSinH;
  }else if(i==7){
    return strElliotSig;
  }else if(i==8){
    return strSoftPlus;
  }else if(i==9){
    return strBentIdentity;
  }else if(i==10){
    return strSinusoid;
  }else if(i==11){
    return strGaussian;
  }else if(i==12){
    return strSinc;
  }else{
    return strRelu;
  }

}

void MakeStrVec(arma::vec nstrVec, String* strVec){
  
  int n = nstrVec.n_elem;
  int nIndex = 1;
  for(int i=1;i<=n;i++){
    nIndex = nstrVec[i-1];
    strVec[i-1] = Num2ActiveStr(nIndex);
  }
  
}





double CrossEntropyVal(arma::vec yVec, arma::vec tVec){
  
  double eps = 1e-7;
  double out = 0;
  //arma::mat tmpMat(1,1);
  //tmpMat = -sum( t %  log(y+eps),0 ); 
  out = -sum( tVec %  log(yVec+eps));
  return out; 
}


arma::mat CrossEntropy(arma::mat y, arma::mat t){
  
  int p = y.n_rows;
  int n = y.n_cols;
  arma::mat out(n,1);
  
  arma::vec yVec(p);
  arma::vec tVec(p);
  
  for(int i=1;i<=n;i++){
    yVec = y.col(i-1);
    tVec = t.col(i-1);
    out(i-1,0) = CrossEntropyVal(yVec, tVec)/n;
  }
  
  return out; 
}

arma::mat NewCrossEntropy(arma::mat y, arma::mat t){
  double eps = 1e-7;
  arma::mat out = -sum( t %  log(y+eps), 0);
  return out.t(); 
}



arma::mat Softmax(arma::mat X){
  
  int p = X.n_rows;
  int n = X.n_cols;
  
  arma::vec x;
  arma::mat out(p,n);
  out.zeros();
  
  double c = 0;
  double sum_exp = 0;
  
  arma::vec exp_x(p);
  
  for(int i = 1;i<=n; i++){
    x = X.col(i-1);
    c = x.max();
    exp_x = exp(x-c);
    sum_exp = sum(exp_x);
    out.col(i-1) = exp_x / sum_exp;
  }
  
  return out;
  
}



arma::mat Masking(arma::mat X, double cut_off){
  
  X.elem( find(X > cut_off) ).ones();
  X.elem( find(X <= cut_off) ).zeros();
  return X;
  
}

arma::mat AlphaMasking(arma::mat X, double cut_off, double alpha){
  
  X.elem( find(X > cut_off) ).ones();
  X.elem( find(X <= cut_off) ).zeros();
  X.elem( find(X <= 0) ) += alpha;
  
  return X;
  
}




arma::vec BatchNorm(arma::vec x){
  
  int n = x.n_elem;
  
  double muB = sum(x)/n;
  
  arma::vec x_mu = x - muB;
  
  double sigB2 = sum( x_mu % x_mu )/n;
  double eps = 1e-7;
  
  double SigEps = sigB2 + eps;
  
  return x_mu/sqrt(SigEps); 
  
}

arma::vec BackwardBatchNorm(arma::vec x, arma::vec dOut){
  
  int n = x.n_elem;
  
  arma::mat GradMat(n,n);
  
  double muB = sum(x)/n;
  
  arma::vec x_mu = x - muB;
  
  double sigB2 = sum( x_mu % x_mu )/n;
  double eps = 1e-7;
  
  double SigEps = sigB2 + eps;
  
  GradMat = ( x_mu*x_mu.t() +  SigEps );
  GradMat.diag() -= n*SigEps;
  
  // for(int i=1;i<=n;i++){
  //   GradMat(i-1,i-1) -= n*SigEps; 
  // }
  // 
  return -1/(n*sqrt(SigEps*SigEps*SigEps)) * GradMat * dOut; 
  
  
}




#endif



















