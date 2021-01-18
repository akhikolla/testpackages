
#ifndef xxBENTIDENTITY_H
#define xxBENTIDENTITY_H


class BentIdentity{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  BentIdentity(){
    n=0;
    p=0;
  }
  
  BentIdentity(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat BentIdentity::Get_Out(){
  return Out;
}

arma::mat BentIdentity::Get_dOut(){
  return dOut;
}

void BentIdentity::forward(arma::mat X){
  
  arma::mat XX = X%X+1;
  Out = ( X + ( sqrt(XX) -1 )/2  ) ;  
  

}


void BentIdentity::backward(arma::mat xX, arma::mat xdOut){
  arma::mat XX = xX%xX+1;
  dOut = xdOut% (1+ 0.5*xX /sqrt(1+xX%xX) ) ;

}


#endif



















