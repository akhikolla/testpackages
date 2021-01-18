
#ifndef xxSIGMOID_H
#define xxSIGMOID_H


class Sigmoid{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  Sigmoid(){
    n=0;
    p=0;
  }
  
  Sigmoid(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xdOut);
  
  
};

arma::mat Sigmoid::Get_Out(){
  return Out;
}

arma::mat Sigmoid::Get_dOut(){
  return dOut;
}

void Sigmoid::forward(arma::mat X){
  
  Out = 1/(1+exp(-X));  
  
  // arma::vec o1(p);
  // o1.ones();
  // 
  // arma::vec x(p);
  // x.zeros();
  // 
  // 
  // for(int i=1;i<=n;i++){
  //   x = X.col(i-1);
  //   Out.col(i-1) = o1/( 1+exp(-x));
  // }

}


void Sigmoid::backward(arma::mat xdOut){

  
  dOut =  ( (1-Out)%Out ) %xdOut;
  // arma::vec x(p);
  // x.zeros();
  // 
  // 
  // 
  // for(int i=1;i<=n;i++){
  //   x = Out.col(i-1);
  //   dOut.col(i-1) =  ( (1-x)%x )  % xdOut.col(i-1) ;
  // }


}


#endif



















