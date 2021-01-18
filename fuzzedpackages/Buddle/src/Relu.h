
#ifndef xxRELU_H
#define xxRELU_H


class Relu{
  
private:
  
  int n;
  int p;
  
  arma::mat Out;
  arma::mat dOut;
  
public:
  
  Relu(){
    n=0;
    p=0;
  }
  
  Relu(int xp, int xn) // Constructor
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

arma::mat Relu::Get_Out(){
  return Out;
}

arma::mat Relu::Get_dOut(){
  return dOut;
}

void Relu::forward(arma::mat X){
  double cut_off=0;
  Out = Masking(X, cut_off) % X;

}


void Relu::backward(arma::mat xdOut){
  
  dOut = Out % xdOut;  
  
  // arma::vec x(p);
  // x.zeros();
  // 
  // for(int i=1;i<=n;i++){
  //   x = Out.col(i-1);
  //   dOut.col(i-1) =  x  % xdOut.col(i-1) ;
  // }
}

#endif



















