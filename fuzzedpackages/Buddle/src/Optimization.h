
#ifndef xxOPTIMIZATION_H
#define xxOPTIMIZATION_H

class Optimization{

private:
  
  int q;
  int p;
  int bRand;
  
  double d_learning_rate;
  double momentum_alpha;
  double d_decay;
  int nIter;
  double d_lr_t;
  double beta1;
  double beta2;
  
  String strOpt;
  
  arma::mat W;
  arma::mat b;
  arma::mat dW;
  arma::mat db;
  
  arma::mat v;
  arma::mat dv;
  
  
  arma::mat vW;
  arma::mat vb;
  arma::mat vv;
  
  arma::mat hW;
  arma::mat hb;
  arma::mat hv;
  
  arma::mat mW;
  arma::mat mb;
  arma::mat mv;
  
  arma::mat nW;
  arma::mat nb;
  arma::mat nv;
  
public:
  
  Optimization(){
    q = 1;
    p = 1;
  }
  
  Optimization(int xq, int xp, int xbRand, double xd_learning_rate, String xstrOpt) // Constructor
    : W(xq,xp), b(xq,1), dW(xq,xp), db(xq,1), v((xq+2),1), dv((xq+2), 1),  
      vW(xq,xp), vb(xq,1), vv((xq+2),1), 
      hW(xq,xp), hb(xq,1), hv((xq+2),1), 
      mW(xq,xp), mb(xq,1), mv((xq+2),1), 
      nW(xq,xp), nb(xq,1), nv((xq+2),1)  { // Default matrix member variable initialization
    
    q = xq;
    p = xp;
    
    bRand = xbRand;
    
    d_learning_rate=xd_learning_rate;
    momentum_alpha = 0.9;
    d_decay = 0.99;
    nIter=0;
    d_lr_t = d_learning_rate;
    beta1=0.9;
    beta2=0.999;
    
    strOpt = xstrOpt;
    
    W.zeros();
    b.zeros();
    dW.zeros();
    db.zeros();
    
    v.zeros();
    dv.zeros();
    
    vW.zeros();
    vb.zeros();
    vv.zeros();
    
    hW.zeros();
    hb.zeros();
    hv.zeros();
    
    mW.zeros();
    mb.zeros();
    mv.zeros();
    
    nW.zeros();
    nb.zeros();
    nv.zeros();
    
    
    

  }
  
  arma::mat Get_W();
  arma::mat Get_b();
  arma::mat Get_v();
  
  
  void Set_W(arma::mat xW);
  void Set_b(arma::mat xb);
  void Set_dW(arma::mat xdW);
  void Set_db(arma::mat xdb);
  
  void Set_v(arma::mat xv);
  void Set_dv(arma::mat xdv);
  
  
  void Update();
  
};

arma::mat Optimization::Get_W(){
  return W;
}


arma::mat Optimization::Get_b(){
  return b;
}

arma::mat Optimization::Get_v(){
  return v;
}


void Optimization::Set_W(arma::mat xW){
  W = xW;
}

void Optimization::Set_b(arma::mat xb){
  b = xb;
}

void Optimization::Set_dW(arma::mat xdW){
  dW = xdW;
}

void Optimization::Set_db(arma::mat xdb){
  db = xdb;
}

void Optimization::Set_v(arma::mat xv){
  v = xv;
}

void Optimization::Set_dv(arma::mat xdv){
  dv = xdv;
}



void Optimization::Update(){
  
  double delta = 1e-5;
  double dtmp = 0;
  
  if(strOpt == strSGD){
    
    W -= d_learning_rate*dW;
    b -= d_learning_rate*db;
    
    if(bRand==1){
      v -= d_learning_rate*dv;
    }
    
    
  }else if(strOpt == strMomentum){
    
    vW = momentum_alpha*vW - d_learning_rate * dW;
    vb = momentum_alpha*vb - d_learning_rate * db;
    
    W += vW;
    b += vb;
    
    if(bRand==1){
      vv = momentum_alpha*vv - d_learning_rate * dv;
      v += vv;
    }
    
    
    
    
  }else if(strOpt == strNesterov){
    
    dtmp = momentum_alpha;
    
    vW = dtmp * vW - d_learning_rate*dW;
    vb = dtmp * vb - d_learning_rate*db;
    
    vW = (dtmp*dtmp)* vW - (1+dtmp)*d_learning_rate * dW;
    vb = (dtmp*dtmp)* vb - (1+dtmp)*d_learning_rate * db;
    
    W += vW;
    b += vb;
    
    if(bRand==1){
      vv = dtmp * vv - d_learning_rate*dv;
      vv = (dtmp*dtmp)* vv - (1+dtmp)*d_learning_rate * dv;
    }
    
    
    
  }else if(strOpt == strAdaGrad){
    
    
    hW += (dW%dW);
    hb += (db%db);
    
    W -= d_learning_rate*(dW/(sqrt(hW)+delta));
    b -= d_learning_rate*(db/(sqrt(hb)+delta));
    
    if(bRand==1){
      hv += (dv%dv);
      v -= d_learning_rate*(dv/(sqrt(hv)+delta));
    }
    
    
    
  }else if(strOpt == strRMSprop){
    
    hW *= d_decay;
    hb *= d_decay;
    
    hW += (1-d_decay)* (dW%dW);
    hb += (1-d_decay)* (db%db);
    
    W -= d_learning_rate*(dW/(sqrt(hW)+delta) );
    b -= d_learning_rate*(db/(sqrt(hb)+delta));
    
    if(bRand==1){
      hv *= d_decay;
      hv += (1-d_decay)* (dv%dv);
      v -= d_learning_rate*(dv/(sqrt(hv)+delta));
    }
    
  }else if(strOpt == strAdam){
    
    nIter++;
    d_lr_t = d_learning_rate*(1-beta2)/(1-beta1);
    
    mW += (1-beta1)*(dW-mW);
    mb += (1-beta1)*(db-mb);
    
    nW += (1-beta2)*(dW%dW - nW);
    nb += (1-beta2)*(db%db - nb);
    
    W -= d_lr_t*mW/sqrt(nW+delta);
    b -= d_lr_t*mb/sqrt(nb+delta);
    
    
    if(bRand==1){
      mv += (1-beta1)*(dv-mv);
      nv += (1-beta2)*(dv%dv - nv);
      v -= d_lr_t*mv/sqrt(nv+delta);
      
    }
    
    
  }else{
    
    W -= d_learning_rate*dW;
    b -= d_learning_rate*db;
    
    if(bRand==1){
      v -= d_learning_rate*dv;
    }
    
  }
  
  v(1,0) = std::abs(v(1,0))+delta;   ////// The second element of v is sigma, and hence should be >0
  
}






#endif



















