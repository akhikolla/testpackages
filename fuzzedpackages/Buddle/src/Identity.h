
#ifndef xxIDENTITY_H
#define xxIDENTITY_H


class Identity{
  
private:
  
  int p;
  int n;
  
  arma::mat Out;
  arma::mat dOut;
   
public:
  
  Identity(){
    n=0;
    p=0;
  }
  
  Identity(int xp, int xn) // Constructor
    : Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    
  }
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat xX);
  void backward(arma::mat xdOut);
  
  
};

arma::mat Identity::Get_Out(){
  return Out;
}

arma::mat Identity::Get_dOut(){
  return dOut;
}

void Identity::forward(arma::mat X){
  
  Out = X;  
  

}


void Identity::backward(arma::mat xdOut){
  dOut = xdOut;

}


#endif



















