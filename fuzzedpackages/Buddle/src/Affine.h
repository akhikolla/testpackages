
#ifndef xxAFFINE_H
#define xxAFFINE_H


class Affine{
  
protected:
  
  int q;     //////  W = q x p matrix,  b = q x 1 vector
  int p;     //////  X = p x n matrix
  int n;     //////  WX + b = q x n matrix, Wx+b = q x 1 vector    
  
  arma::mat Out;   //// q x n
  arma::mat dOut;  ////  _dOut = the same dimension with Out, q x n
                   ////  dOut = the same dimension with X,  p x n
  
private:
  
  arma::mat dW;    //// the same dimension with W, q x p
  arma::mat db;    //// the same dimension with b, q x 1
  
  arma::mat W;
  arma::mat b;

  int bRand;
  FInv finv;
  String strDist;
  
  arma::mat v;
  arma::mat dv;
  
  
public:
  
  Affine(){
    q=1;
    p=1;
    n=1;
  }
  
  Affine(int xq, int xp, int xn, int xbRand, String xstrDist) // Constructor
    : Out(xq, xn), dOut(xp, xn), dW(xq, xp), db(xq, 1),
      W(xq, xp), b(xq,1), finv(FInv(xq, 1, xstrDist)), v((xq+2),1), dv((xq+2),1) { // Default matrix member variable initialization
    
    
    q = xq;
    p = xp;
    n = xn;
    
    bRand = xbRand;
    strDist = xstrDist;
    
    Out.zeros();
    dOut.zeros();
    
    W.zeros();
    b.zeros();
    
    dW.zeros();
    db.zeros();
    
    v.zeros();
    dv.zeros();
    
  }
  
  
  
  arma::mat Get_Out();
  arma::mat Get_dOut();
  
  arma::mat Get_dW();
  arma::mat Get_db();
  
  arma::mat Get_dv();
  
  
  void Set_W(arma::mat xW);
  void Set_b(arma::mat xb);
  
  void Set_v(arma::mat xv);
  
  
  void forward(arma::mat xX);
  void backward(arma::mat xX, arma::mat xdOut);
  
  
};

arma::mat Affine::Get_Out(){
  return Out;
}

arma::mat Affine::Get_dOut(){
  return dOut;
}

arma::mat Affine::Get_dW(){
  return dW;
}

arma::mat Affine::Get_db(){
  return db;
}

arma::mat Affine::Get_dv(){
  return dv;
}


void Affine::Set_W(arma::mat xW){
  W = xW;
}

void Affine::Set_b(arma::mat xb){
  b = xb;
}

void Affine::Set_v(arma::mat xv){
  v = xv;
  
  v.rows(0,1) = xv.rows(0,1);
  
  arma::vec x(q);
  x.randu(q);
  
  for(int i=1;i<=q;i++){
    v(i-1+2,0) = x[i-1];
  }
  
}


void Affine::forward(arma::mat xX){
  
  arma::mat OneMat(1, n);
  OneMat.ones();

  Out = W*xX+ b*OneMat;
  
  if(bRand==1){
    
    finv.forward(v);
    Out += finv.Get_Out()*OneMat;
  }
  

}


void Affine::backward(arma::mat xX, arma::mat xdOut){

  dW = xdOut* xX.t();
  db = sum(xdOut, 1);
  dOut = W.t()* xdOut;
  
  if(bRand==1){
    
    finv.backward(v, db);   //// db or sum(_dOut, 1);
    dv = finv.Get_dOut();
    
  }
  
  

}

class gAffine: public Affine {
  
private:
  
  int q2;
  arma::mat V;
  arma::mat dV;
  
  arma::mat tmp_Out;
  arma::mat tmp_dOut;
  
  Link link;
  FInv finv2;
  
  String strLink;
  
public:
  
  gAffine(){
    q=1;
    p=1;
    n=1;
  }
  
  gAffine( int xq, int xp, int xn, int xbRand, String xstrDist, String xstrLink) // Constructor
    : Affine(xq, xp, xn, xbRand, xstrDist), V((xq+2), xn), dV((xq+2), xn), 
      tmp_Out(xq,n), tmp_dOut( (xq+2),n), link(Link(xq, xn, xstrLink)), finv2(FInv(xq, n, xstrLink)){
      
      strLink = xstrLink;
      
      q2 = xq+2;  
      V.zeros();
      dV.zeros();
      
      tmp_Out.zeros();
      tmp_dOut.zeros();
  }
  
  arma::mat Get_dV();
  void Set_V(arma::mat xV);
  
  void gforward(arma::mat xX);
  void gbackward(arma::mat xX, arma::mat xdOut);
  

};

void gAffine::Set_V(arma::mat xV){
  
  V.rows(0,1) = xV.rows(0,1);
  
  arma::mat tmpU(q,n);
  tmpU.randu(q,n);
  
  V.rows(3-1, q2-1) = tmpU;
  
}



void gAffine::gforward(arma::mat xX){
  
  forward(xX);            ///// Get Out  after affine layer
  tmp_Out = Out;          ////// Out:qxn                      
  link.forward(Out);
  V = link.Get_Out();     ///// Pass Out to link function to make V: q2xn
  finv2.forward(V);        ///// Pass Out to FInv node 
  Out = finv2.Get_Out();      ///// Out:qxn
  
  
}


void gAffine::gbackward(arma::mat xX, arma::mat xdOut){
                                 ////  _dOut :qxn 
  finv2.backward(V, xdOut);    ////   Get dOut after FInv
  tmp_dOut = finv2.Get_dOut();   //// tmp_dOut: q2xn
  link.backward(tmp_Out, tmp_dOut);
  
  xdOut = link.Get_dOut();          /////_dOut:qxn
  
  backward(xX, xdOut);        ///// Pass _dOut to affine layer 
                              /////  Get dW, db, dV
  
}



#endif



















