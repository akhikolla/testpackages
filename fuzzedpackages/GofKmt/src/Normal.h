
#ifndef xxNORMAL_H
#define xxNORMAL_H


class Normal{
  
private:
  
  int n;
  arma::vec X;
  
  arma::vec SubGiVec1;
  arma::vec SubGiVec2;
  arma::vec SubGiVec3;
  
  double Div;
  int PM;
  int PM2;
  
public:
  
  Normal(){
    Div=1;
    PM=1;
    PM2=1;
    
  }
  
  Normal(arma::vec xX, arma::mat IntMat) 
      : X(xX.n_elem), SubGiVec1(2561), SubGiVec2(2561), SubGiVec3(2561){
    
    n = xX.n_elem;
    X = xX;
    
    Div=128;
    PM=1280;
    PM2=2561;
    
    
    for(int i=1;i<=3;i++){
      if(i==1){
        SubGiVec1 = IntMat.col(i-1);
      }else if(i==2){
        SubGiVec2 = IntMat.col(i-1);
      }else{
        SubGiVec3 = IntMat.col(i-1);
      }
    }

    
    
  }
  
  
  double fn(double x);
  double Fn(double x);
  
  double Gam11(double x);
  double Gam12(double x);
  double Gam13(double x);
  
  double Gam23(double x);
  double Gam22(double x);
  double Gam33(double x);
  
  arma::mat GammaMatrix(double x);
  
  double phix(double x);
  
  double c1(double x);
  double c2(double x);
  double c3(double x);
  
  double subgi1(double y);
  double subgi2(double y);
  double subgi3(double y);
  double gi(double y, double Xi);
  
  double subGi(double x, int nI);
  double Gi(double y, double Xi);
  
  
  
  
};





double Normal::fn(double x){
  return R::dnorm(x, 0, 1, 0);
}

double Normal::Fn(double x){
  return R::pnorm(x, 0, 1, 1, 0);
}


double Normal::phix(double x){
  return x;
}


double Normal::Gam11(double x){
  double out= 1-Fn(x);
  return out;
}


double Normal::Gam12(double x){
  double out= fn(x);
  return out;
}


double Normal::Gam13(double x){
  double out= x*fn(x);
  return out;
}


double Normal::Gam23(double x){
  double out= (x*x+1)*fn(x);
  return out;
}

double Normal::Gam22(double x){
  
  double out=0;
  out = x*fn(x)+(1-Fn(x));
  return out;
}

double Normal::Gam33(double x){
  
  double out=0;
  out = (x*x*x+x)*fn(x)+2*(1-Fn(x));
  return out;
  
}

arma::mat Normal::GammaMatrix(double x){
  
  arma::mat out(3,3);
  
  double fnVal=fn(x);
  double FnVal=Fn(x);
  double xfnVal = x*fnVal;
  out(0,0) = 1-FnVal; //Gam11(x);
  out(0,1) = fnVal; //Gam12(x);
  out(0,2) = xfnVal; //Gam13(x);
  
  out(1,0) = fnVal;
  out(1,1) = x*fnVal+(1-FnVal);            // Gam22(x);
  out(1,2) = (x*x+1)*fnVal; //Gam23(x);
  
  out(2,0) = xfnVal;
  out(2,1) = out(1,2);
  out(2,2) = (x*x*x+x)*fnVal+2*(1-FnVal);    //Gam33(x);
  
  return out;
}


double Normal::c1(double x){
  double out=0;
  double fnVal = fn(x);
  double FnVal = 1 - Fn(x);
  out = -(x*x+1)*fnVal*fnVal + (x*x*x+3*x)*fnVal*FnVal + 2*FnVal*FnVal;
  
  return out;
  
}

double Normal::c2(double x){
  
  double out=0;
  double fnVal = fn(x);
  double FnVal = 1 - Fn(x);
  out = 2*FnVal*FnVal*FnVal + (x*x*x+3*x)*fnVal*FnVal*FnVal - (2*x*x+3)*fnVal*fnVal*FnVal + x*fnVal*fnVal*fnVal;
  
  return out;
  
} 


double Normal::c3(double x){
  double out=0;
  double fnVal = fn(x);
  double FnVal = 1 - Fn(x);
  
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  
  out = (4-x2)*fnVal*fnVal -(2*x3+10*x)* fnVal*FnVal  - (x4+3)*FnVal*FnVal;
  
  return out;
  
}



double Normal::subgi1(double y){
  
  double out=0;
  double fnVal = fn(y);
  double FnVal = 1 - Fn(y);
  out = 2*fnVal*( FnVal*FnVal+y*fnVal*FnVal - fnVal*fnVal  )/ c2(y);
  
  return out;
}


double Normal::subgi2(double y){
  
  double out=0;
  double fnVal = fn(y);
  double fnVal2 = fnVal*fnVal;
  double fnVal3 = fnVal2*fnVal;
  double fnVal4 = fnVal3*fnVal;
  
  double FnVal = 1 - Fn(y);
  double FnVal2 = FnVal*FnVal;
  double FnVal3 = FnVal2*FnVal;
  double FnVal4 = FnVal3*FnVal;
  
  double y2 = y*y;
  double y3 = y2*y;
  double y4 = y3*y;
  double y5 = y4*y;
  
  out = fnVal*( 4*y*FnVal4 + (2*y4+8*y2-2)*fnVal*FnVal3 +
    (y5-7*y)*fnVal2*FnVal2 - (2*y4+3*y2-1)*fnVal3*FnVal +
    (y3+y)*fnVal4 )/ (c1(y)*c2(y));
  
  return out;
}



double Normal::subgi3(double y){
  
  double out=0;
  double fnVal = fn(y);
  double fnVal2 = fnVal*fnVal;
  double fnVal3 = fnVal2*fnVal;
  double fnVal4 = fnVal3*fnVal;
  
  double FnVal = 1 - Fn(y);
  double FnVal2 = FnVal*FnVal;
  double FnVal3 = FnVal2*FnVal;
  double FnVal4 = FnVal3*FnVal;
  
  double y2 = y*y;
  double y3 = y2*y;
  double y4 = y3*y;
  double y5 = y4*y;
  
  out = fnVal*( 2*(y2-1)*FnVal4 + (y5+2*y3-9*y)*fnVal*FnVal3 -
    (4*y4+9*y2-5)*fnVal2*FnVal2 + (5*y3+9*y)*fnVal3*FnVal -
    2*(y2+1)*fnVal4 )/ (c1(y)*c2(y));
  
  return out;
}

double Normal::gi(double y, double Xi){
  
  double out=0;
  out = subgi1(y) + Xi*subgi2(y) + (Xi*Xi-1)*subgi3(y);
  
  return out;    
}



double Normal::subGi(double x, int nI){
  
  
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
  
  //int nLen = 100000;
  
  double out = 0;
  double xi = SP;
  
  //arma::mat tmpMat(3,1);
  
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


double Normal::Gi(double y, double Xi){
  
  double out=0;
  
  out = subGi(y,1) + Xi*subGi(y,2) + (Xi*Xi-1)*subGi(y,3);
  
  
  return out;    
}




#endif



















