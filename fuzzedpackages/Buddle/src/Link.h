
#ifndef xxLINK_H
#define xxLINK_H


class Link{
  
private:
  
  int q;     //////  v: (q+2)x1 or (q+1) vectpr
  int n;
  int q2;
  
  
  String strLink;
  
  arma::mat Out;   //// q x 1
  arma::mat dOut;  ////  _dOut = the same dimension with Out, q x n
                   ////  dOut = the same dimension with X,  p x n
  
public:
  
  Link(){
    q=1;
    n=1;
  }
  
  Link(int xq, int xn, String xstrLink ) // Constructor
    : Out((xq+2), xn), dOut(xq, xn) { // Default matrix member variable initialization
    
    q = xq;
    n = xn;
    q2 = q+2;
    
    strLink = xstrLink;

    Out.zeros();
    dOut.zeros();

  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();

  void forward(arma::mat X);
  void backward(arma::mat X, arma::mat xdOut);
  
  
};

arma::mat Link::Get_Out(){
  return Out;
}

arma::mat Link::Get_dOut(){
  return dOut;
}



void Link::forward(arma::mat X){
  
  double del = 1e-5;                     /// X:qxn, Out:q2xn, dOut: qxn
  double mu,sig;
  arma::vec x(q);
  arma::mat tmpU(q,n);
  tmpU.randu(q,n);
  double pi = 3.14159;
  for(int i=1;i<=n;i++){
    
    x = X.col(i-1);                /// x:qx1
    mu = sum(x)/q;
    
    if(strLink == strNormal){
      sig = sqrt( sum((x-mu)%(x-mu)) /q )+del;
      Out(0, i-1) = mu;
      Out(1, i-1) = sig;
      
    }else if(strLink == strLogistic){
      sig = sqrt(3/pi)* sqrt( sum((x-mu)%(x-mu)) /q )+del;
      Out(0, i-1) = mu;
      Out(1, i-1) = sig;
      
    }else if(strLink == strGamma){
      sig = mu;
      Out(0, i-1) = mu;
      Out(1, i-1) = sig;
      
    }else{
      
      sig = sqrt( sum((x-mu)%(x-mu)) /q )+del;
      Out(0, i-1) = mu;
      Out(1, i-1) = sig;
      
    }
    
  }
  
  Out.rows(3-1, q2-1) = tmpU;
  
}

void Link::backward(arma::mat X, arma::mat xdOut){

  double del=1e-5;                        /// X:qxn   xdOut: q2xn
  double mu,sig;                          //// dOut: qxn, Out:q2xn
  arma::vec z(q);
  double d1,d2;
  
  for(int i=1;i<=n;i++){
    mu = Out(0, i-1);
    sig = Out(1, i-1)+del;
    
    if(strLink == strNormal){
      z = ( X.col(i-1)-mu )/sig;
      d1 = xdOut(0, i-1);
      d2 = xdOut(1, i-1);
      dOut.col(i-1) = (d1+d2*z)/q;
    }else if(strLink == strPoisson){
      z.ones();
      d1 = xdOut(0, i-1);
      d2 = xdOut(1, i-1);
      dOut.col(i-1) = (d1+d2*z/(2*sig))/q;
      
    }else if(strLink == strGamma){
      z.ones();
      d1 = xdOut(0, i-1);
      d2 = xdOut(1, i-1);
      dOut.col(i-1) = (d1+d2*z)/q;
      
    }else{
      z = ( X.col(i-1)-mu )/sig;
      d1 = xdOut(0, i-1);
      d2 = xdOut(1, i-1);
      dOut.col(i-1) = (d1+d2*z)/q;
    }

  }

}


#endif



















