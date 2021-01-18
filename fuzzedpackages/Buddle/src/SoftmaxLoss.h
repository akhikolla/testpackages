
#ifndef xxSOFTMAXLOSS_H
#define xxSOFTMAXLOSS_H

class SoftmaxLoss{

private:
  
  int r;     //////  y = r x n matrix,  t = r x n matrix
  int n;     //////  loss = n x 1 matrix
  
  double loss;
  
  arma::mat Entropy;   ////  n x 1
  arma::mat y;  
  
  arma::mat dOut;   //// dOut = r x n matrix
                   
  
public:
  
  SoftmaxLoss(){
    r = 1;
    n = 1;
  }
  
  SoftmaxLoss(int xr, int xn) // Constructor
    : Entropy(xn,1), y(xr, xn), dOut(xr,xn) { // Default matrix member variable initialization
    
    r = xr;
    n = xn;
    
    loss= 0;
    
    Entropy.zeros();
    y.zeros();
    dOut.zeros();
  }
  
  arma::mat Get_Entropy();
  arma::mat Get_y();
  arma::mat Get_dOut();
  double Get_loss();
  
  void forward(arma::mat X, arma::mat xt); 
  void forward_predict(arma::mat X); 
  void backward(arma::mat xt);
  
};

arma::mat SoftmaxLoss::Get_Entropy(){
  return Entropy;
}


arma::mat SoftmaxLoss::Get_y(){
  return y;
}


arma::mat SoftmaxLoss::Get_dOut(){
  return dOut;
}

double SoftmaxLoss::Get_loss(){
  return loss;
}

void SoftmaxLoss::forward(arma::mat X, arma::mat xt){
  
  double eps = 1e-7;
  y = Softmax(X);
  arma::mat out = -sum( xt %  log(y+eps), 0);
  Entropy = out.t();
  loss = accu(Entropy)/n;
}

void SoftmaxLoss::forward_predict(arma::mat X){
  
  y = Softmax(X);
  
}





void SoftmaxLoss::backward(arma::mat xt){

  dOut = (y - xt) / n;
  
}





#endif



















