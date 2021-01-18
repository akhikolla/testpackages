
#ifndef xxLAYER_H
#define xxLAYER_H


class Layer{
  
protected:  
  
  int q;
  int p;
  int n;
  int bAct;
  int bBatch;
  int bDrop;
  int bTest;
  double drop_ratio;
  double d_learning_rate;
  double d_initial_weight;
  
  int bRand;
  String strDist;
  
  
  String strAct;          ///  Activation function
  String strOpt;          ///   Optimization method 
  
  Affine affine;
  Relu relu;
  Sigmoid sigmoid;
  LeakyRelu leakyrelu;
  TanH tanh;
  ArcTan arctan;
  ArcSinH arcsinh;
  ElliotSig elliotsig;
  SoftPlus softplus;
  BentIdentity bentidentity;
  Sinusoid sinusoid;
  Identity identity;
  Gaussian gaussian;
  Sinc sinc;
  
  
  Dropout dropout;
  Batchnorm batchnorm;
  
  arma::mat W;
  arma::mat b;
  arma::mat dW;
  arma::mat db;
  
  arma::mat v;
  arma::mat dv;
  
  Optimization opt;
  arma::mat Out;
  arma::mat dOut;
  

public:
  
  Layer(){
    q=1;
    p=1;
    n=1;
  }
  
  Layer(int xq, int xp, int xn, int xbAct, int xbBatch, int xbDrop, int xbTest, double xdrop_ratio, 
        double xd_learning_rate, double xd_initial_weight, String xstrAct, String xstrOpt, int xbRand, String xstrDist)
    : affine(Affine(xq, xp, xn, xbRand, xstrDist)), relu(Relu(xq,xn)), sigmoid(Sigmoid(xq,xn)), leakyrelu(LeakyRelu(xq,xn)),
      tanh(TanH(xq,xn)), arctan(ArcTan(xq,xn)), arcsinh(ArcSinH(xq,xn)),
      elliotsig(ElliotSig(xq,xn)), softplus(SoftPlus(xq,xn)),
      bentidentity(BentIdentity(xq,xn)), sinusoid(Sinusoid(xq,xn)), identity(Identity(xq,xn)),
      gaussian(Gaussian(xq,xn)), sinc(Sinc(xq,xn)),
      dropout( Dropout(xq, xn, xbTest, xdrop_ratio)), batchnorm(Batchnorm(xq,xn)),
      W(xq,xp), b(xq,1), dW(xq,xp), db(xq,1), v((xq+2),1), dv((xq+2),1),
      opt(Optimization(xq, xp, xbRand, xd_learning_rate, xstrOpt)),
      Out(xq,xn), dOut(xp,xn) { // Default matrix member variable initialization
    
    
    bRand=xbRand;
    strDist = xstrDist;
    
    q = xq;
    p = xp;
    n = xn;
    
    
    bAct = xbAct;
    bBatch=xbBatch;
    bDrop=xbDrop;
    bTest=xbTest;
    drop_ratio=xdrop_ratio;
    d_learning_rate = xd_learning_rate;
    d_initial_weight = xd_initial_weight;
    
    strOpt = xstrOpt;
    strAct = xstrAct; 
     
    W.randn(q, p);
    b.zeros();
    
    v.randu((q+2),1);
    v(1,0) += 1e-5;
    
    W *= d_initial_weight;
    b *= d_initial_weight;
    double dp2 = p/2;
    double dp = p;
    
    if(strAct == strRelu || strAct == strLeakyRelu){
      dp2 = std::sqrt(dp2);
      W /= dp2;
    }else{
      dp = std::sqrt(dp);
      W /= dp;
    }
    
    
    
  }
  
  
  arma::mat Get_Out();
  arma::mat Get_dOut();  
  
  arma::mat Get_W();
  arma::mat Get_b();  
  
  arma::mat Get_dW();
  arma::mat Get_db();  
  
  arma::mat Get_v();
  arma::mat Get_dv();  
  
  
  void Set_W(arma::mat xW);
  void Set_b(arma::mat xb);
  
  void Set_v(arma::mat xv);
  
  
  void forward(arma::mat X);
  void backward(arma::mat X, arma::mat dOut);
  
  
};

arma::mat Layer::Get_Out(){
  return Out;
}

arma::mat Layer::Get_dOut(){
  return dOut;
}

arma::mat Layer::Get_W(){
  
  return W;
}


arma::mat Layer::Get_b(){
  return b;
}

arma::mat Layer::Get_dW(){
  return dW;
}

arma::mat Layer::Get_db(){
  return db;
}


arma::mat Layer::Get_v(){
  return v;
}

arma::mat Layer::Get_dv(){
  return dv;
}




void Layer::Set_W(arma::mat xW){
  
  W = xW;
}


void Layer::Set_b(arma::mat xb){
  
  b=xb;
}

void Layer::Set_v(arma::mat xv){
  
  v=xv;
}


void Layer::forward(arma::mat X){
  

  affine.Set_W(W);
  affine.Set_b(b);
  affine.Set_v(v);
  
  affine.forward(X);
  
  Out = affine.Get_Out();
  
  if(bAct == 1){
    
    if(bBatch==1){
      batchnorm.forward(Out);
      Out = batchnorm.Get_Out();
    }
    
    if(strAct == strRelu){
      relu.forward(Out);
      Out = relu.Get_Out();
    }else if(strAct == strSigmoid){
      sigmoid.forward(Out);
      Out = sigmoid.Get_Out();
    }else if(strAct == strLeakyRelu){
      leakyrelu.forward(Out);
      Out = leakyrelu.Get_Out();
    }else if(strAct == strTanH){
      tanh.forward(Out);
      Out = tanh.Get_Out();
      
    }else if(strAct == strArcTan){
      arctan.forward(Out);
      Out = arctan.Get_Out();
      
    }else if(strAct == strArcSinH){
      arcsinh.forward(Out);
      Out = arcsinh.Get_Out();
      
    }else if(strAct == strElliotSig){
      elliotsig.forward(Out);
      Out = elliotsig.Get_Out();
      
    }else if(strAct == strSoftPlus){
      softplus.forward(Out);
      Out = softplus.Get_Out();
      
    }else if(strAct == strBentIdentity){
      bentidentity.forward(Out);
      Out = bentidentity.Get_Out();
      
    }else if(strAct == strSinusoid){
      sinusoid.forward(Out);
      Out = sinusoid.Get_Out();
      
    }else if(strAct == strGaussian){
      gaussian.forward(Out);
      Out = gaussian.Get_Out();
      
    }else if(strAct == strSinc){
      sinc.forward(Out);
      Out = sinc.Get_Out();
      
    }else{
      identity.forward(Out);
      Out = identity.Get_Out();
    }
    
    if(bDrop==1){
      dropout.forward(Out);
      Out = dropout.Get_Out();
    }
    
    
  }
  
}


void Layer::backward(arma::mat xX, arma::mat xdOut){

  if(bAct == 1){
    if(bDrop==1){
      dropout.backward(xdOut);
      xdOut = dropout.Get_dOut();
    }
    
    if(strAct == strRelu){
      relu.backward(xdOut);
      xdOut = relu.Get_dOut();
    }else if(strAct == strSigmoid){
      sigmoid.backward(xdOut);
      xdOut = sigmoid.Get_dOut();
    }else if(strAct == strLeakyRelu){
      leakyrelu.backward(xdOut);
      xdOut = leakyrelu.Get_dOut();
    }else if(strAct == strTanH){
      tanh.backward(xdOut);
      xdOut = tanh.Get_dOut();
    }else if(strAct == strArcTan){
      
      if(bDrop==1){
        arctan.backward(dropout.Get_Out(), xdOut);
      }else{
        arctan.backward(affine.Get_Out(), xdOut);
      }
        xdOut = arctan.Get_dOut();
    }else if(strAct == strArcSinH){
      
      if(bDrop==1){
        arcsinh.backward(dropout.Get_Out(), xdOut);
      }else{
        arcsinh.backward(affine.Get_Out(), xdOut);
      }
      xdOut = arcsinh.Get_dOut();
    }else if(strAct == strElliotSig){
      
      if(bDrop==1){
        elliotsig.backward(dropout.Get_Out(), xdOut);
      }else{
        elliotsig.backward(affine.Get_Out(), xdOut);
      }
      xdOut = elliotsig.Get_dOut();
    }else if(strAct == strSoftPlus){
      
      if(bDrop==1){
        softplus.backward(dropout.Get_Out(), xdOut);
      }else{
        softplus.backward(affine.Get_Out(), xdOut);
      }
      xdOut = softplus.Get_dOut();
    }else if(strAct == strBentIdentity){
      
      if(bDrop==1){
        bentidentity.backward(dropout.Get_Out(), xdOut);
      }else{
        bentidentity.backward(affine.Get_Out(), xdOut);
      }
      xdOut = bentidentity.Get_dOut();
    }else if(strAct == strSinusoid){
      
      if(bDrop==1){
        sinusoid.backward(dropout.Get_Out(), xdOut);
      }else{
        sinusoid.backward(affine.Get_Out(), xdOut);
      }
      xdOut = sinusoid.Get_dOut();
    }else if(strAct == strGaussian){
      
      if(bDrop==1){
        gaussian.backward(dropout.Get_Out(), xdOut);
      }else{
        gaussian.backward(affine.Get_Out(), xdOut);
      }
      xdOut = gaussian.Get_dOut();
    }else if(strAct == strSinc){
      
      if(bDrop==1){
        sinc.backward(dropout.Get_Out(), xdOut);
      }else{
        sinc.backward(affine.Get_Out(), xdOut);
      }
      xdOut = sinc.Get_dOut();
    }else{
      identity.backward(xdOut);
      xdOut = identity.Get_dOut();
      
    }
    
    
    if(bBatch==1){
      batchnorm.backward(xdOut);
      xdOut = batchnorm.Get_dOut();
    }
    
  }
  
  affine.backward(xX, xdOut);
  dOut = affine.Get_dOut();
  
  dW = affine.Get_dW();
  db = affine.Get_db();
  dv = affine.Get_dv();  
  ////////////// Optimization
  
  opt.Set_W(W);
  opt.Set_b(b);
  opt.Set_dW(dW);
  opt.Set_db(db);
  opt.Update();
  
  W = opt.Get_W();
  b = opt.Get_b();
  
  
}





#endif



















