
#ifndef xxMOMENTUM_H
#define xxMOMENTUM_H


class Momentum{
  
private:
  
  double d_learning_rate;
  
  
public:
  
  Momentum(){
    d_learning_rate=0;
  }
  
  Momentum(double xd_learning_rate) { // Default matrix member variable initialization
    
    d_learning_rate = xd_learning_rate;
    
  }

  arma::mat Update(arma::mat W, arma::mat dW);
  
};


arma::mat Momentum::Update(arma::mat W, arma::mat dW, arma::mat vMat){
  
  int q = W.n_rows;
  int p = W.n_cols;
  arma::mat Out(q,p);
  
  vMat = Momentum_alpha*vMat - d_learning_rate * dW;
  Out = W+vMat;
  
  return Out;
  
}



#endif



















