
#ifndef xxTANH_H
#define xxTANH_H


class TanH{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  TanH(){
    n=0;
    p=0;
  }
  
  TanH(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xdOut);
  
  
};

arma::mat TanH::Get_Out(){
  return Out;
}

arma::mat TanH::Get_dOut(){
  return dOut;
}

void TanH::forward(arma::mat X){
  
  Out = (exp(X)-exp(-X))/(exp(X)+exp(-X));  
  

}


void TanH::backward(arma::mat xdOut){

  
  dOut =  (1 -  Out%Out ) % xdOut;


}


#endif



















