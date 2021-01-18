
#ifndef xxSOFTPLUS_H
#define xxSOFTPLUS_H


class SoftPlus{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  SoftPlus(){
    n=0;
    p=0;
  }
  
  SoftPlus(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat SoftPlus::Get_Out(){
  return Out;
}

arma::mat SoftPlus::Get_dOut(){
  return dOut;
}

void SoftPlus::forward(arma::mat X){
  
  Out = log(1+exp(X)) ;
  

}


void SoftPlus::backward(arma::mat xX, arma::mat xdOut){
  dOut = xdOut/(1+exp(- xX));

}


#endif



















