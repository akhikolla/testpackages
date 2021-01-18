
#ifndef xxDROPOUT_H
#define xxDROPOUT_H


class Dropout{
  
private:
  
  int n;
  int p;
  int bTest;
  double drop_ratio;
  
  arma::mat Mask;
  arma::mat Out;
  arma::mat dOut;
  
public:
  
  Dropout(){
    n=0;
    p=0;
    bTest=0;
    drop_ratio = 0;
  }
  
  Dropout(int xp, int xn, int xbTest, double xdrop_ratio) // Constructor
    : Mask(xp, xn), Out(xp, xn), dOut(xp, xn) { // Default matrix member variable initialization
    
    n = xn;
    p = xp;
    bTest = xbTest;
    drop_ratio = xdrop_ratio;
    
    Mask.zeros();
    Out.zeros();
    dOut.zeros();
    
  }
  
  arma::mat Get_Mask();
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  void forward(arma::mat X);
  void backward(arma::mat xdOut);
  
  
};

arma::mat Dropout::Get_Mask(){
  return Mask;
}


arma::mat Dropout::Get_Out(){
  return Out;
}

arma::mat Dropout::Get_dOut(){
  return dOut;
}

void Dropout::forward(arma::mat X){
  
  arma::mat tmpMat(p,n);
  tmpMat.randu();
  
  if(bTest==0){
    Mask = Masking(tmpMat, drop_ratio);
    Out = X % Mask;
  }else{
    Out = (1-drop_ratio)*X; 
  }
  
}


void Dropout::backward(arma::mat xdOut){
  
  dOut = Mask % xdOut;  
  
}

#endif



















