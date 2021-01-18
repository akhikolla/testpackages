
#ifndef xxLOGISTIC_H
#define xxLOGISTIC_H


class Logistic{
  
private:
  
  int n;
  arma::vec X;
  
  arma::vec SubGiVec1;
  arma::vec SubGiVec2;
  arma::vec SubGiVec3;
  
  double Div;
  int PM;
  int PM2;
  
  arma::vec ReVec;
  double ReDiv;
  int RePM;
  int RePM2;
  
  
  double Oldx;
  int rstar;
  int rstar2;
  
  
public:
  
  Logistic(){
    Div=1;
    PM=1;
    PM2=1;
    
    Oldx=-20;
    rstar=1;
    rstar2=1;
    
    
  }
  
  Logistic(arma::vec xX, arma::mat IntMat, arma::mat ReMat) 
    : X(xX.n_elem), SubGiVec1(1281), SubGiVec2(1281), SubGiVec3(1281), ReVec(1281){
    
    n = xX.n_elem;
    X = xX;
    
    Div=64;
    PM=640;
    PM2=1281;
    
    for(int i=1;i<=3;i++){
      if(i==1){
        SubGiVec1 = IntMat.col(i-1);
      }else if(i==2){
        SubGiVec2 = IntMat.col(i-1);
      }else{
        SubGiVec3 = IntMat.col(i-1);
      }
    }
    
    
    //////// Re part
    ReDiv=64;
    RePM=640;
    RePM2=1281;
    
    for(int j=1;j<=(RePM2);j++){
      
      ReVec[j-1] = ReMat(j-1,0);
    }
    
    
    Oldx=-20;
    rstar=1;
    rstar2=1;
    
    
    
    
  }
  
  
  double fn(double x);
  double Fn(double x);
  double phix(double x);
  
  double re(double x);
  double Re(double x);
  arma::mat GammaMatrix(double x);
  
  double g1(double x);
  double g2(double x);
  
  double k1(double x);
  double k2(double x);
  
  double cofx(double x);
  double A12A22A21(double x);
  double B11(double x);
  
  
  
  double subgi1(double y);
  double subgi2(double y);
  double subgi3(double y);
  double gi(double y, double Xi);
  
  double subGi(double x, int nI);
  double Gi(double y, double Xi);
  
  void Set_rstar();
  
};


void Logistic::Set_rstar(){
  rstar=1;
  rstar2=1;
  Oldx=-20;
}

double Logistic::fn(double x){
  return R::dlogis(x, 0, 1, 0);
}

double Logistic::Fn(double x){
  return R::plogis(x, 0, 1, 0, 0);
}

double Logistic::phix(double x){
  double ex = exp(x);
  return ((ex-1)/(ex+1));
}

double Logistic::re(double x){
  double ex = exp(x);
  double x2 = x*x;
  double exm1 = (1-ex);
  double exp1 = (1+ex);
  
  return x2*ex*exm1*exm1/(exp1*exp1*exp1*exp1);
}

double Logistic::Re(double x){
  
  int nIndex = 0;
  double SP = 0;
  double dInc = 1/(ReDiv); 
  
  
  if(x>Oldx){
    rstar=rstar2;
  }else if(x<Oldx){
    rstar=1;
  }
  
  
  
  for(int ith=rstar;ith <= (RePM2 - 1); ith++){
    if(x< -15){
      nIndex = -1;
      SP = -16;
      break;
    }
    
    if(x >= 15){
      
      nIndex=(RePM2 - 1);
      SP=15;
      break;
    }
    if( (x >= -15 + (ith-1)*dInc )&(x < -15 + ith*dInc ) ){
      rstar2 = ith;
      
      nIndex  = ith-1;
      SP = -15 + (ith-1)*dInc;
      break;
    }
    
  }
  
  
  Oldx = x;
  
  if(nIndex == -1){
    return ReVec[0];
    
  }
  
  if(nIndex == (RePM2 - 1)){
    return ReVec[nIndex];
  }
  
  const double nGap = 1e-3;  // 10^5
  
  int nLen = 1e+3;     

  double out = 0;
  double xi = SP;

  double dResid = 0;
  
  for(int i=1; i<=nLen; i++){
    xi += nGap;
    
    if(xi>=x){
      
      dResid = re(x) *(nGap + x - xi) ;
      break;
    }
  
    out += re(xi);
    
  }
  
  return ReVec[nIndex] - (out*nGap + dResid);
  
  
}


arma::mat Logistic::GammaMatrix(double x){
  
  arma::mat out(3,3);
  
  double ex = exp(x);
  double e2x = pow(ex, 2);
  double ex1 = ex+1;
  double ex1_3 = pow(ex1, 3);
  
  double flVal = fn(x);
  double FlVal = Fn(x);
  
  double logVal = log(1+ex)/3;
  double Rex = Re(x);
  
  out(0,0) = 1-FlVal;
  out(0,1) = flVal;
  out(0,2) = x*flVal;
  out(1,0) = out(0,1);
  out(1,1) = (3*e2x + 1)/(3*ex1_3);
  out(1,2) =  logVal - flVal* (x*(3+e2x) + ex1)/(3*ex1);
  out(2,0) = out(0,2);
  out(2,1) = out(1,2);
  out(2,2) = FlVal - 2*x*flVal + Rex-1;
  return out;
}



double Logistic::g1(double x){
  double ex = exp(x);
  double e2x = ex*ex;
  
  double flVal = fn(x);
  double logVal = log(1+ex)/3;
  
  return logVal - flVal* (x*(3+e2x) + 1+ex)/(3*(1+ex));
}


double Logistic::g2(double x){
  
  double Rex = Re(x);
  
  double flVal = fn(x);
  double FlVal = Fn(x);
  
  return FlVal - 2*x*flVal + Rex-1;
  
}

double Logistic::k1(double x){
  
  double ex = exp(x);
  double g1Val = g1(x);
  double g2Val = g2(x);
  
  double ex1 = (1+ex);
  
  return 3*ex*ex1*ex1*ex1*g2Val - 3*x*ex*ex1*ex1*ex1*g1Val;
  
}

double Logistic::k2(double x){
  
  double ex = exp(x);
  double e2x = ex*ex;
  double g1Val = g1(x);
  
  double ex1 = (1+ex);
  return -3*ex*ex1*ex1*ex1*g1Val+x*ex*(1+3*e2x);
  
}

double Logistic::cofx(double x){
  
  double ex = exp(x);
  double e2x = ex*ex;
  double g1Val = g1(x);
  double g2Val = g2(x);
  
  double ex1 = (1+ex);
  double ex1_3 = ex1*ex1*ex1;
  double ex1_6 = ex1_3*ex1_3;
  
  return 1/(3*(1+3*e2x)*ex1_3*g2Val - 9* ex1_6*g1Val*g1Val);
}


double Logistic::A12A22A21(double x){
  
  double cofxVal = cofx(x);
  
  double ex = exp(x);
  double e2x = ex*ex;
  double g1Val = g1(x);
  double g2Val = g2(x);
  
  double ex1 = (1+ex);
  double ex1_2 = ex1*ex1;
  double ex1_3 = ex1*ex1*ex1;
  
  return (9*cofxVal*e2x*ex1_2* ( 3*ex1_3*g2Val- 6*x*ex1_3*g1Val +x*x*(1+3*e2x)) );
}

double Logistic::B11(double x){
  
  double AAA = A12A22A21(x);
  double ex = exp(x);
  double ex1 = (1+ex);
  double ex1_2 = ex1*ex1;
  
  return 1/(3*ex1_2 - AAA);
}


double Logistic::subgi1(double x){
  
  double cofxVal = cofx(x);
  
  double ex = exp(x);
  double B11Val = B11(x);
  double k1Val = k1(x);
  double k2Val = k2(x);
  
  return 3*ex*(1+ex)*B11Val* (1+3*cofxVal*( (1-ex)*k1Val + (1+ex +x-x*ex )*k2Val ) );
  
  
}

double Logistic::subgi2(double x){
  
  double cofxVal = cofx(x);
  
  double ex = exp(x);
  //double e2x = ex*ex;
  double g1Val = g1(x);
  double g2Val = g2(x);
  
  double ex1 = (1+ex);
  double ex1_2 = ex1*ex1;
  double ex1_3 = ex1*ex1*ex1;
  
  double B11Val = B11(x);
  double k1Val = k1(x);
  double k2Val = k2(x);
  
  double ans =  (-9*ex*ex1_2*B11Val*cofxVal*k1Val) - 3*ex*cofxVal*(3*ex1_3*(1-ex)*g2Val - 3*ex1_3*g1Val*(ex1+x-x*ex) ) -
    B11Val*27*ex*ex1_2*cofxVal*cofxVal*((1-ex)*k1Val*k1Val+ (ex1+x-x*ex)*k1Val*k2Val  );
  return ans;
}


double Logistic::subgi3(double x){
  
  double cofxVal = cofx(x);
  
  double ex = exp(x);
  double e2x = ex*ex;
  double g1Val = g1(x);
  //double g2Val = g2(x);
  
  double ex1 = (1+ex);
  double ex1_2 = ex1*ex1;
  double ex1_3 = ex1*ex1*ex1;
  
  double B11Val = B11(x);
  double k1Val = k1(x);
  double k2Val = k2(x);
  
  double ans1 = (  -9*ex*ex1_2* B11Val*cofxVal * k2Val );
  double ans2 = -3*ex*cofxVal*( -3*ex1_3*(1-ex )*g1Val + (1+3*e2x)* (1+ex+x-x*ex) );
  double ans3 = - B11Val*27*ex* ex1_2* cofxVal*cofxVal * ( (1-ex)*k1Val*k2Val +(1+ex+x-x*ex)* (k2Val*k2Val)   );
  
  return ans1 + ans2+ans3;
  
}



double Logistic::gi(double y, double Xi){
  
  double out=0;
  out = subgi1(y) + phix(Xi)*subgi2(y) + (Xi*phix(Xi)-1)*subgi3(y);
  
  return out;    
}




double Logistic::subGi(double x, int nI){
  
  int nIndex = 0;
  double SP = 0;
  double dInc = 1/(Div); 
  
  for(int ith=1;ith <= (PM2 - 1); ith++){
    if(x<-10){
      nIndex = -1;
      SP = -11;
      break;
    }
    
    if(x >= 10){
      
      nIndex=(PM2 - 1);
      SP=10;
      break;
    }
    if( (x >= -10 + (ith-1)*dInc )&(x < -10 + ith*dInc ) ){
      nIndex  = ith-1;
      SP = -10 + (ith-1)*dInc;
      break;
    }
    
  }

  if(nIndex == -1){
    
    if(nI == 1){
      return SubGiVec1[0];
    }else if(nI == 2){
      return SubGiVec2[0];
    }else{
      return SubGiVec3[0];
    }
    
  }
  
  if(nIndex == (PM2 - 1)){
    if(nI == 1){
      return SubGiVec1[nIndex];
    }else if(nI == 2){
      return SubGiVec2[nIndex];
    }else{
      return SubGiVec3[nIndex];
    }
  }
  
  const double nGap = 1e-4;  // 10^5
  
  int nLen = 1e+5;     

  double out = 0;
  double xi = SP;

  double dResid = 0;
  
  for(int i=1; i<=nLen; i++){
    xi += nGap;
    
    if(xi>=x){
      
      if(nI==1){
        dResid = subgi1(x) *(nGap + x - xi) ;
      }else if(nI==2){
        dResid = subgi2(x) *(nGap + x - xi) ;
      }else{
        dResid = subgi3(x) *(nGap + x - xi) ;
      }
      
      break;
    }
    
    
    if(nI == 1){
      out += subgi1(xi);
      
    }else if(nI == 2){
      out += subgi2(xi);
      
    }else{
      out += subgi3(xi);
      
    }
    
    
  }
  
  if(nI == 1){
    return out*nGap + SubGiVec1[nIndex]+dResid;
  }else if(nI == 2){
    return out*nGap + SubGiVec2[nIndex] +dResid ;
  }else{
    return out*nGap + SubGiVec3[nIndex] +dResid;
  }
  
}



double Logistic::Gi(double y, double Xi){
  
  double out=0;
  
  out = subGi(y,1) + phix(Xi)*subGi(y,2) + (Xi*phix(Xi)-1)*subGi(y,3);
  
  
  return out;    
}






#endif



















