
#ifndef xxL2LOSS_H
#define xxL2LOSS_H

class L2loss{

private:
  
  int r;     //////  y = r x n matrix,  t = r x n matrix
  int n;     //////  loss = n x 1 matrix
  
  double loss;
  
  arma::mat y;  
  arma::mat dOut;   //// dOut = r x n matrix

public:
  
  L2loss(){
    r = 1;
    n = 1;
  }
  
  L2loss(int xr, int xn) // Constructor
    : y(xr, xn), dOut(xr,xn) { // Default matrix member variable initialization
    
    r = xr;
    n = xn;
    loss = 0;
    
    y.zeros();
    dOut.zeros();
  }
  
  arma::mat Get_y();
  arma::mat Get_dOut();
  double Get_loss();
  
  void forward(arma::mat X, arma::mat xt); 
  void forward_predict(arma::mat X); 
  void backward(arma::mat xt);
  
};



arma::mat L2loss::Get_y(){
  return y;
}


arma::mat L2loss::Get_dOut(){
  return dOut;
}

double L2loss::Get_loss(){
  
  return loss;
}


void L2loss::forward(arma::mat X, arma::mat xt){
  
  y = X;
  
  arma::mat yt = abs(y - xt);
  loss = 0.5* accu( yt%yt)/n;
  
}

void L2loss::forward_predict(arma::mat X){
  
  y = X;
  
}


void L2loss::backward(arma::mat xt){

  dOut = (y - xt)/n;
  
}



#endif



















