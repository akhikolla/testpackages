
#ifndef xxELLIOTSIG_H
#define xxELLIOTSIG_H


class ElliotSig{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  ElliotSig(){
    n=0;
    p=0;
  }
  
  ElliotSig(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat ElliotSig::Get_Out(){
  return Out;
}

arma::mat ElliotSig::Get_dOut(){
  return dOut;
}

void ElliotSig::forward(arma::mat X){
  
  Out = X / (1+abs(X));
  

}


void ElliotSig::backward(arma::mat xX, arma::mat xdOut){
  dOut = xdOut/  ((1+ abs(xX)) % (1+ abs(xX))) ;

}


#endif



















