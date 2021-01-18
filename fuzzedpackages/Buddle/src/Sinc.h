
#ifndef xxSINC_H
#define xxSINC_H


class Sinc{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  Sinc(){
    n=0;
    p=0;
  }
  
  Sinc(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat Sinc::Get_Out(){
  return Out;
}

arma::mat Sinc::Get_dOut(){
  return dOut;
}

void Sinc::forward(arma::mat X){
  
  double del = 1e-7;
  arma::uvec ind = find(X == 0);
  
  X.ones();
  Out = X%sin(X) / (X+del);
  Out.elem( ind ).ones();

}


void Sinc::backward(arma::mat X, arma::mat xdOut){
  
  double del = 1e-7;
  arma::uvec ind = find(X == 0);
  
  dOut = xdOut % ( cos(X)/(X+del) - sin(X)/( X%X+del ));
  dOut.elem( ind ).zeros();
  

}


#endif



















