
#ifndef xxCAUCHY_H
#define xxCAUCHY_H


class Cauchy{
  
private:
  
  int n;
  arma::vec X;
  
  arma::vec SubGiVec1;
  arma::vec SubGiVec2;
  arma::vec SubGiVec3;
  
  double Div;
  int PM;
  int PM2;
  
  double pi;
  
  
public:
  
  Cauchy(){
    Div=1;
    PM=1;
    PM2=1;
    
    
  }
  
  Cauchy(arma::vec xX, arma::mat IntMat) 
    : X(xX.n_elem), SubGiVec1(2561), SubGiVec2(2561), SubGiVec3(2561){
    
    n = xX.n_elem;
    X = xX;
    
    Div=128;
    PM=1280;
    PM2=2561;
    
    
    pi = 3.14159265;
    
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
  double phix(double x);
  
  ////// Gamma part
  double Gam11(double x);
  double Gam12(double x);
  double Gam13(double x);
  double Gam22(double x);
  double Gam23(double x);
  double Gam33(double x);
  arma::mat GammaMatrix(double x);
  
  double Dc(double x);
  double Bc(double x);
  double Ac(double x);
  
  
  double subgi1(double y);
  double subgi2(double y);
  double subgi3(double y);
  double gi(double y, double Xi);
  
  double subGi(double x, int nI);
  double Gi(double y, double Xi);
  
  
  
};



double Cauchy::fn(double x){
  return ( 1 / ( pi*(1+pow(x,2)) ));
}

double Cauchy::Fn(double x){
  return R::pcauchy(x,0,1,1,0);
}

double Cauchy::phix(double x){
  
  return ((2*x)/(1+x*x));
}



double Cauchy::Gam11(double x){
  
  double x21 = 1+pow(x,2);
  double denom = 4*pi*pow(x21,2);
  double out = 2*pi*pow(x21,2) - 4*pow(x21,2)*atan(x);
  return out/denom;
}

double Cauchy::Gam12(double x){
  
  double x21 = 1+pow(x,2);
  double denom = 4*pi*pow(x21,2);
  double out= 4*x21;
  return out/denom;
}

double Cauchy::Gam13(double x){
  
  double x21 = 1+pow(x,2);
  double denom = 4*pi*pow(x21,2);
  double out= 4*x*x21;
  return out/denom;
}

double Cauchy::Gam23(double x){
  double x21 = 1+pow(x,2);
  double denom = 4*pi*pow(x21,2);
  
  double out= 4*pow(x,2);
  return out/denom;
}


double Cauchy::Gam22(double x){
  
  double x21 = 1+pow(x,2);
  double denom = 4*pi*pow(x21,2);
  double x2_1 = 1-pow(x,2);
  
  double out=0;
  out = (pi-2*atan(x))*pow(x21,2)+2*x*x2_1;
  return out/denom;
}


double Cauchy::Gam33(double x){
  
  double x21 = 1+pow(x,2);
  double denom = 4*pi*pow(x21,2);
  double x2_1 = 1-pow(x,2);
  
  double out=0;
  out = (pi-2*atan(x))*pow(x21,2)-2*x*x2_1;
  return out/denom;
  
}


arma::mat Cauchy::GammaMatrix(double x){
  
  arma::mat out(3,3);
  
  out(0,0) = Gam11(x);
  out(0,1) = Gam12(x);
  out(0,2) = Gam13(x);
  
  out(1,0) = out(0,1);
  out(1,1) = Gam22(x);
  out(1,2) = Gam23(x);
  
  out(2,0) = out(0,2);;
  out(2,1) = out(1,2);
  out(2,2) = Gam33(x);
  
  return out;
}


double Cauchy::Dc(double x){
  
  double out= (pi-2*atan(x));
  
  return out;
  
}


double Cauchy::Bc(double x){
  
  double x21 = 1+pow(x,2);
  double dc = Dc(x);
  
  double out=x21*dc - 2*x;
  
  return out;
  
} 


double Cauchy::Ac(double x){
  
  double x21 = 1+pow(x,2);
  double dc = Dc(x);
  
  double out=x21*pow(dc,2) +2*x*dc-8;
  
  return out;
  
}



double Cauchy::subgi1(double x){
  
  double x21 = 1+pow(x,2);
  double bc = Bc(x);
  double ac = Ac(x);
  double denom = x21*ac*bc;
  
  
  double out=2*pow(bc,2)/denom;
  
  return out;
}


double Cauchy::subgi2(double x){
  
  double x21 = 1+pow(x,2);
  double bc = Bc(x);
  double dc = Dc(x);
  double ac = Ac(x);
  double denom = x21*ac*bc;
  
  
  double out= (-8*bc+8*x*(x21*pow(dc, 2)-4))/denom;
  
  return out;
}

double Cauchy::subgi3(double x){
  
  double out=0;
  double x21 = 1+pow(x,2);
  double bc = Bc(x);
  double dc = Dc(x);
  double ac = Ac(x);
  double denom = x21*ac*bc;
  
  
  out = ( -8*x*bc+4*(pow(x,4)-1)*pow(dc,2) - 8*x*x21*dc+32)/denom;  
  
  return out;
}


double Cauchy::gi(double y, double Xi){
  
  double out=0;
  out = subgi1(y) + phix(Xi)*subgi2(y) + (Xi*phix(Xi)-1)*subgi3(y);
  
  return out;    
}


double Cauchy::subGi(double x, int nI){
  
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



double Cauchy::Gi(double y, double Xi){
  
  double out=0;
  
  out = subGi(y,1) + phix(Xi)*subGi(y,2) + (Xi*phix(Xi)-1)*subGi(y,3);
  
  
  return out;    
}






#endif



















