
#ifndef xxSINUSOID_H
#define xxSINUSOID_H


class Sinusoid{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  Sinusoid(){
    n=0;
    p=0;
  }
  
  Sinusoid(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat Sinusoid::Get_Out(){
  return Out;
}

arma::mat Sinusoid::Get_dOut(){
  return dOut;
}

void Sinusoid::forward(arma::mat X){
  
  Out = sin(X);  
  

}


void Sinusoid::backward(arma::mat xX, arma::mat xdOut){
  dOut = xdOut% cos(xX) ;

}


#endif



















