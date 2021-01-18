
#ifndef xxGAUSSIAN_H
#define xxGAUSSIAN_H


class Gaussian{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  Gaussian(){
    n=0;
    p=0;
  }
  
  Gaussian(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat Gaussian::Get_Out(){
  return Out;
}

arma::mat Gaussian::Get_dOut(){
  return dOut;
}

void Gaussian::forward(arma::mat X){
  
  arma::mat XX = -X%X;
  
  Out = exp(XX);  
}


void Gaussian::backward(arma::mat xX, arma::mat xdOut){
  
  arma::mat XX = -xX%xX;
  dOut = 2*xdOut % xX % exp(XX);

}


#endif



















