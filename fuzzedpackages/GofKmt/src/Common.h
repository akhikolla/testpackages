
#ifndef xxCOMMON_H
#define xxCOMMON_H

//#include <RcppArmadillo.h>


String strNormal("Normal");
String strLogistic("Logistic");
String strCauchy("Cauchy");
String strExponential("Exponential");

String strPoisson("Poisson");
String strGamma("Gamma");




double AbsVal(double x);
arma::vec GetLineVec(int n);


arma::vec GetLineVec(arma::vec X, int nNum){
  
  
  double dAdd = 2.5;
  int n=X.n_elem;
  int nLen = (n+1)*nNum;
  
  arma::vec out(nLen);
  
  int nIndex=0;
  
  double SP=0;
  double EP=0;
  double del = 1e-3;
  
  double dGap = 0;
  for(int i=1;i<=n;i++){
    
    nIndex = (i-1)*nNum;
    
    if(i==1){
      EP = X[i-1]-del;
      SP = EP - dAdd;
    }else{
      EP = X[i-1]-del;
      SP = X[i-2];
    }
    dGap = (EP-SP)/(nNum-1);
    
    for(int j=1;j<=nNum;j++){
      out[nIndex+j-1] = SP + dGap*(j-1);
    }
    
  }
  
  
  ////////////// After Xn
  nIndex = n*nNum;
  
  SP = X[n-1];
  EP = SP + dAdd;
  
  dGap = (EP-SP)/(nNum-1);
  
  for(int j=1;j<=nNum;j++){
    out[nIndex+j-1] = SP + dGap*(j-1);
  }
  
  
  
  return out;
}

double AbsVal(double x){
  
  if(x<0){
    return -x;
  }else{
    return x;
  }
  
}





#endif



















