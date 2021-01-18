#include "Dist2.h"
//#include "armadillo"


Dist2::Dist2(int _n, int _m, int _nModel2){

  n = _n;
  m = _m;
  nm = n*m;
  nModel2 = _nModel2;

  psi = new double[n];
}

Dist2::~Dist2(){
  delete [] psi;
}


void Dist2::Setpsi(){

  ///// Normal Case

  if(nModel2 == 1){
    for(int i=1;i<=n;i++){
      psi[i-1]=0;
    }
  }else{
    for(int i=1;i<=n;i++){
      psi[i-1]=1;
    }
  }

  return;

}


double* Dist2::Getpsi(){
  return psi;
}



