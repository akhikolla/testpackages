
#ifndef xxLEAKYRELU_H
#define xxLEAKYRELU_H


class LeakyRelu{
  
private:
  
  int n;
  int p;
  
  arma::mat Out;
  arma::mat dOut;
  
public:
  
  LeakyRelu(){
    n=0;
    p=0;
  }
  
  LeakyRelu(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;

    Out.zeros();
    dOut.zeros();
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat X);
  void backward(arma::mat xdOut);
  
  
};

arma::mat LeakyRelu::Get_Out(){
  return Out;
}

arma::mat LeakyRelu::Get_dOut(){
  return dOut;
}

void LeakyRelu::forward(arma::mat X){
  double cut_off=0;
  double alpha = 0.1;
  Out = AlphaMasking(X, cut_off, alpha) % X;

}


void LeakyRelu::backward(arma::mat xdOut){
  
  dOut = Out % xdOut;  
  
  
}

#endif



















