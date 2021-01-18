
#ifndef xxARCSINH_H
#define xxARCSINH_H


class ArcSinH{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  ArcSinH(){
    n=0;
    p=0;
  }
  
  ArcSinH(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat ArcSinH::Get_Out(){
  return Out;
}

arma::mat ArcSinH::Get_dOut(){
  return dOut;
}

void ArcSinH::forward(arma::mat X){
  
  arma::mat XX = 1+X%X;
  Out = log( X + sqrt( XX )  );  
  

}


void ArcSinH::backward(arma::mat xX, arma::mat xdOut){
  
  arma::mat XX = 1+xX%xX;
  dOut = xdOut/ sqrt(XX) ;

}


#endif



















