
#ifndef xxARCTAN_H
#define xxARCTAN_H


class ArcTan{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  ArcTan(){
    n=0;
    p=0;
  }
  
  ArcTan(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat ArcTan::Get_Out(){
  return Out;
}

arma::mat ArcTan::Get_dOut(){
  return dOut;
}

void ArcTan::forward(arma::mat X){
  
  Out = atan(X);  
  

}


void ArcTan::backward(arma::mat xX, arma::mat xdOut){
  dOut = xdOut/(1+xX%xX) ;

}


#endif



















