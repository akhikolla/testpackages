
#ifndef xxKMT_H
#define xxKMT_H



class Kmt{
  
private:
  
  int n;
  arma::vec X;
  arma::mat GiMat;
  arma::vec T2;
  arma::vec Modified_T2;
  
  String strDist;
  
  Normal norm1;
  Logistic logis1;
  Cauchy cauchy1;
  
  double OptimalX;
  double OptimalFVal;
  
  double OptimalXP;
  double OptimalFValP;
  double OptimalXM;
  double OptimalFValM;
  


public:
  
  Kmt(){
    n = 1 ; 
    
  }
  
  Kmt(arma::vec xX, int xn, arma::mat NormalMat, arma::mat LogisMat, arma::mat ReMat, arma::mat CauchyMat, String xstrDist)
   : X(xn), GiMat(xn, xn), T2(3), norm1(Normal(xX, NormalMat)), logis1(Logistic(xX, LogisMat, ReMat)), 
     cauchy1(Cauchy(xX, CauchyMat)){
    
    GiMat.zeros();
    T2.zeros();
    n = xn;
    X = xX;
    
    strDist = xstrDist;
    
    
    OptimalX = 0;
    OptimalFVal=0;
    
    OptimalXP = 0;
    OptimalFValP=0;
    
    OptimalXM = 0;
    OptimalFValM=0;
    
    
  }

  arma::mat GetGiMat();
  double GetOptX();
  double GetOptFVal();
  arma::vec GetT2();
  
  ////////////////////////
  double GetOptXP();
  double GetOptFValP();
  
  double GetOptXM();
  double GetOptFValM();
  ////////////////////////
  
  
  
  arma::mat GetGammaMatrix(double x);
  void SetGiMat();
  void SetT2();
  
  double hiz(double z, int i);
  double SecantZero(int i, int Mos);
  
  double ObjVal(double z);
  void FindOptimal();
  
  
  double RawObjVal(double z);
  void Modified_SetT2();
  void Modified_FindOptimal();
  
};



arma::mat Kmt::GetGiMat(){
  return GiMat;
}

double Kmt::GetOptX(){
  return OptimalX;
}

double Kmt::GetOptFVal(){
  return OptimalFVal;
}

arma::vec Kmt::GetT2(){
  return T2;
}



double Kmt::GetOptXP(){
  return OptimalXP;
}

double Kmt::GetOptFValP(){
  return OptimalFValP;
}

double Kmt::GetOptXM(){
  return OptimalXM;
}

double Kmt::GetOptFValM(){
  return OptimalFValM;
}






arma::mat Kmt::GetGammaMatrix(double x){
  
  arma::mat out(3,3);
  
  if(strDist==strNormal){
    out = norm1.GammaMatrix(x);
  }else if(strDist == strLogistic){
    out = logis1.GammaMatrix(x);
  }else{
    out = cauchy1.GammaMatrix(x);
  }
  
  return out;
}


void Kmt::SetGiMat(){
  
  arma::mat prod1(n,3);
  arma::mat prod2(3,n);
  
  double Xi=0;
  
  for(int i=1;i<=n;i++){
    Xi = X[(i-1)];
    prod1(i-1,0)=1;
    
    if(strDist == strNormal){
      
      prod1(i-1,1)=norm1.phix(Xi);
      prod1(i-1,2)=(Xi*norm1.phix(Xi)-1);
      
      prod2(0,i-1) = norm1.subGi(Xi,1);
      prod2(1,i-1) = norm1.subGi(Xi,2);
      prod2(2,i-1) = norm1.subGi(Xi,3);
      
      
    }else if(strDist == strLogistic){
      
      prod1(i-1,1)=logis1.phix(Xi);
      prod1(i-1,2)=(Xi*logis1.phix(Xi)-1);
      
      
      prod2(0,i-1) = logis1.subGi(Xi,1);
      prod2(1,i-1) = logis1.subGi(Xi,2);
      prod2(2,i-1) = logis1.subGi(Xi,3);
      
      
      
    }else{
      
      prod1(i-1,1)=cauchy1.phix(Xi);
      prod1(i-1,2)=(Xi*cauchy1.phix(Xi)-1);
      
      
      prod2(0,i-1) = cauchy1.subGi(Xi,1);
      prod2(1,i-1) = cauchy1.subGi(Xi,2);
      prod2(2,i-1) = cauchy1.subGi(Xi,3);
      
      
    }
    
    
  }
  GiMat = prod1*prod2;
 
  /////////////////////////////////
  logis1.Set_rstar();
  /////////////////////////////////
 
 
}


double Kmt::hiz(double z, int i){
  
  double out=0;
  double Xi = 0;
  
  for(int k=(i+1);k<=n;k++){
    Xi = X[k-1];
    
    if(strDist == strNormal){
      out -= norm1.gi(z, Xi);
    }else if(strDist == strLogistic){
      out -= logis1.gi(z, Xi);
    }else{
      out -= cauchy1.gi(z, Xi);
    }
    
    
  }
  
  return out;
  
}



double Kmt::SecantZero(int i, int Mos){
  
  int nLen = 10000;
  
  double SP=0;
  double EP=0;
  
  if(Mos == 1){
    SP = X[(i-1)];
    EP = X[i];
  }else if(Mos == 2){
    SP = (X[(i-1)] + X[i])/2;
    EP = X[i];
  }else{
    SP = X[i];
    EP = X[(i-1)];
  }
  
  double nGap = (EP-SP)/nLen;
  
  double crit = 0.001;
  int nIter = 500;
  
  double x0 = SP;
  double x1 = x0 + nGap;
  
  double xn1 = 0;
  double xn2 = 0;
  double xn = 0;
  
  double fVal = 0;
  double fVal1 = 0;
  double fVal2 = 0;
  
  double OptimalZero = 0;
  double del = 1e-5;
  
  for(int iter=1;iter<=nIter;iter++){
    
    if(iter == 1){
      xn2 = x0;
      xn1 = x1;
    }
    
    fVal1 = hiz(xn1, i);
    fVal2 = hiz(xn2, i);
    
    if(fVal1 == fVal2){
      xn = xn1 - fVal1*(xn1-xn2)/(fVal1-fVal2+del);
    }else{
      xn = xn1 - fVal1*(xn1-xn2)/(fVal1-fVal2);
    }
    
    fVal = hiz(xn, i);
    
    if(AbsVal(fVal) < crit){
      OptimalZero = xn;
      break;
    }
    xn2 = xn1;
    xn1 = xn;
  }
  
  
  
  
  return OptimalZero;
}


double Kmt::ObjVal(double z){
  
  double out = 0;

  int nIndex = 0;
  
  for(int i=1;i<=n;i++){
    if(z < X[(i-1)]){
      nIndex = i-1;
      break;
    }
  }
  
  if(X[(n-1)] <= z){
    nIndex = n;
  }
  
  double Xi = 0;
  if(nIndex == 0){
    
    for(int i=1;i<=n;i++){
      Xi = X[(i-1)];
      
      if(strDist == strNormal){
        out -= norm1.Gi(z, Xi);
      }else if(strDist == strLogistic){
        out -= logis1.Gi(z, Xi);
      }else{
        out -= cauchy1.Gi(z, Xi);
      }
      
      
    }
    
  }else if(nIndex == n){
    
    for(int i=1;i<=n;i++){
      Xi = X[(i-1)];
      out -= GiMat(i-1,i-1); //Gi(z, Xi)
    }
    out += n;
  }else{
    for(int i=(nIndex+1);i<=n;i++){
      Xi = X[(i-1)];
      
      if(strDist == strNormal){
        out -= norm1.Gi(z, Xi);
      }else if(strDist == strLogistic){
        out -= logis1.Gi(z, Xi);
      }else{
        out -= cauchy1.Gi(z, Xi);
      }
      
    }
    
    for(int i=1;i<=nIndex;i++){
      Xi = X[(i-1)];
      out += (1 - GiMat(i-1,i-1));
    }
  }
  
  double dn = n;
  
  double dOut = 0;
  dOut = AbsVal(out)/ (std::sqrt(dn)) ;
  
  return dOut;
}



double Kmt::RawObjVal(double z){
  
  double out = 0;
  
  int nIndex = 0;
  
  for(int i=1;i<=n;i++){
    if(z < X[(i-1)]){
      nIndex = i-1;
      break;
    }
  }
  
  if(X[(n-1)] <= z){
    nIndex = n;
  }
  
  double Xi = 0;
  if(nIndex == 0){
    
    for(int i=1;i<=n;i++){
      Xi = X[(i-1)];
      
      if(strDist == strNormal){
        out -= norm1.Gi(z, Xi);
      }else if(strDist == strLogistic){
        out -= logis1.Gi(z, Xi);
      }else{
        out -= cauchy1.Gi(z, Xi);
      }
      
      
    }
    
  }else if(nIndex == n){
    
    for(int i=1;i<=n;i++){
      Xi = X[(i-1)];
      out -= GiMat(i-1,i-1); //Gi(z, Xi)
    }
    out += n;
  }else{
    for(int i=(nIndex+1);i<=n;i++){
      Xi = X[(i-1)];
      
      if(strDist == strNormal){
        out -= norm1.Gi(z, Xi);
      }else if(strDist == strLogistic){
        out -= logis1.Gi(z, Xi);
      }else{
        out -= cauchy1.Gi(z, Xi);
      }
      
    }
    
    for(int i=1;i<=nIndex;i++){
      Xi = X[(i-1)];
      out += (1 - GiMat(i-1,i-1));
    }
  }
  
  double dn = n;
  
  double dOut = 0;
  dOut = out/ (std::sqrt(dn)) ;
  
  return dOut;
}




void Kmt::SetT2(){


  double U1=0;
  double U2=0;
  double U=0;
  
  
  double dMax = 0;
  double Optimal = X[0];
  double Xi = 0;
  int nStatus = 0;
  
  
  
  
  for(int i=1;i<=n;i++){
    
    Xi =  X[i-1]; 
    
    U1 = 0;
    
    for(int k=1;k<=n;k++){
      if(k<=i){
        U1 -= GiMat(k-1,k-1);
      }else{
        U1 -= GiMat(k-1,i-1);
      }  
    }
    
    U1 += i; 
    U2 = U1-1;
    
    if(AbsVal(U1) >= AbsVal(U2)){
      U = AbsVal(U1);
      nStatus = 1;
    }else{
      U = AbsVal(U2);
      nStatus = 0;
    }
    
    if(U>dMax){
      dMax = U;
      Optimal = Xi;
    }
    
  
  }
  
  double dn = n;
  
  T2[0] = Optimal;
  T2[1] = nStatus;
  T2[2] = dMax/ ( std::sqrt(dn) );
  
  
 
}




void Kmt::FindOptimal(){

  double Xi=0;
  double Xi1=0;
  double FVal=0;
  
  double OldFVal = T2[2];
  OptimalFVal=OldFVal;
  OptimalX = T2[0];
  double Zerox = 0;
  
  double T2X = 0;
  double T2Y = 0;
  
  //////////////////////
  
  /////////////////////// i=0
  
  Xi = X[0];
  
  Zerox = SecantZero(1, 3);
  
  if( (Zerox < Xi)  | (Zerox >= Xi1) ){
    Zerox = Xi;
  }
  
  if(Zerox != Xi){
    FVal = ObjVal(Zerox);
    
    if(FVal > OldFVal){
      OptimalX = Zerox;
      OptimalFVal = FVal;
      OldFVal = FVal;
    }
    
    
  }
  
  double tmpVal=0;
  double tmpVal2=0;
  double nGap = 0;

  for(int i=1;i<=(n-1);i++){
    Xi = X[i-1];
    Xi1 = X[i];
    
    nGap = (Xi1-Xi)/100;

    tmpVal = hiz( (Xi+nGap), i);
    tmpVal2 = hiz( (Xi1-nGap), i);

    if( (tmpVal*tmpVal2)<0 ){
      
      for(int nMos=1;nMos<=3;nMos++){
        
        Zerox = SecantZero(i, nMos);
        
        if( (Zerox < Xi)  | (Zerox >= Xi1) ){
          Zerox = Xi;
        }
        
        if(Zerox != Xi){
          FVal = ObjVal(Zerox);
          
          if(FVal>T2Y){
            T2Y = FVal;
            T2X = Zerox;
          }
          
          if(FVal > OptimalFVal){
            
            OptimalX = T2X;
            OptimalFVal = FVal;
            
          }
        }
        
      }
      
      
    }
      

  }
    
    
  
  
  
  
  // i=n
  FVal = ObjVal(X[n-1]);
  
  if(FVal>OldFVal){
    OptimalX = X[n-1];
    OptimalFVal = FVal;
    OldFVal = FVal;
  }
  
  
}




void Kmt::Modified_SetT2(){
  
  
  double U1=0;
  double U2=0;
  
  double PU=0;
  double MU=0;
  double OptimalP = X[0];
  double OptimalM = X[0];
  
  double Xi = 0;
  
  
  
  for(int i=1;i<=n;i++){
    
    Xi =  X[i-1]; 
    
    U1 = 0;
    
    for(int k=1;k<=n;k++){
      if(k<=i){
        U1 -= GiMat(k-1,k-1);
      }else{
        U1 -= GiMat(k-1,i-1);
      }  
    }
    
    U1 += i; 
    U2 = U1-1;
    
    if(U1<0){
      if(U1<MU){
        MU=U1;
        OptimalM = Xi;
      }
    }else{
      if(U1>PU){
        PU=U1;
        OptimalP = Xi;
      }
    }
    
    
    
    if(U2<0){
      if(U2<MU){
        MU=U2;
        OptimalM = Xi;
      }
    }else{
      if(U2>PU){
        PU=U2;
        OptimalP = Xi;
      }
    }
    
  }
  
  double dn = n;
  
  OptimalXP = OptimalP;
  OptimalFValP=PU/ ( std::sqrt(dn) );
  
  OptimalXM = OptimalM;
  OptimalFValM=MU/ ( std::sqrt(dn) );
  
  
  
}




void Kmt::Modified_FindOptimal(){
  
  double Xi=0;
  double Xi1=0;
  double FVal=0;
  
  double Zerox = 0;
  
  
  //////////////////////
  
  /////////////////////// i=0
  
  Xi = X[0];
  
  Zerox = SecantZero(1, 3);
  
  if( (Zerox < Xi)  | (Zerox >= Xi1) ){
    Zerox = Xi;
  }
  
  if(Zerox != Xi){
    FVal = RawObjVal(Zerox);
    
    if(FVal > OptimalFValP){
      OptimalXP = Zerox;
      OptimalFValP = FVal;
      
    }
    
    if(FVal < OptimalFValM){
      OptimalXM = Zerox;
      OptimalFValM = FVal;
      
    }
    
    
  }
  
  double tmpVal=0;
  double tmpVal2=0;
  double nGap = 0;
  
  for(int i=1;i<=(n-1);i++){
    Xi = X[i-1];
    Xi1 = X[i];
    
    nGap = (Xi1-Xi)/100;
    
    tmpVal = hiz( (Xi+nGap), i);
    tmpVal2 = hiz( (Xi1-nGap), i);
    
    //////////////////////////
    
    
    /////////////////////////
    
    if( (tmpVal*tmpVal2)<0 ){
      
      
      for(int nMos=1;nMos<=3;nMos++){
        
        Zerox = SecantZero(i, nMos);
        
        if( (Zerox < Xi)  | (Zerox >= Xi1) ){
          Zerox = Xi;
        }
        
        if(Zerox != Xi){
          FVal = RawObjVal(Zerox);
          
          if(FVal > OptimalFValP){
            OptimalXP = Zerox;
            OptimalFValP = FVal;
            
          }
          
          if(FVal < OptimalFValM){
            OptimalXM = Zerox;
            OptimalFValM = FVal;
            
          }
        }
        
      }
      
    }
    
    
  }
  
  
  
  
  // i=n
  FVal = RawObjVal(X[n-1]);
  
  if(FVal > OptimalFValP){
    OptimalXP = X[n-1];
    OptimalFValP = FVal;
    
  }
  
  if(FVal < OptimalFValM){
    OptimalXM = X[n-1];
    OptimalFValM = FVal;
    
  }
  
  
}




#endif



















