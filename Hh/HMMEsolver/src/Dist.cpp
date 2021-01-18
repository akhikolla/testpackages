#include <stdio.h>
#include <math.h>
#include "Dist.h"
//#include "armadillo"


Dist::Dist(int _n, int _m, int _nModel){

  n = _n;
  m = _m;
  nm = n*m;
  nModel = _nModel;

  W = new double[nm];
  s = new double[nm];
  eta = new double[nm];
  deta_dmu = new double[nm];
}

Dist::~Dist(){
  delete [] W;
  delete [] s;
  delete [] eta;
  delete [] deta_dmu;
}


void Dist::SetW(double* _mu){

  ///// Normal Case
  int nInd=0;

  for(int i=1;i<=n;i++){
    for(int j=1;j<=m;j++){
      nInd = (i-1)*m+(j-1);

      if(nModel == 1){
        W[nInd] = 1;
      }else{
        W[nInd] = _mu[nInd];
      }

    }

  }

  return;
}


void Dist::Seteta(double* _mu){

  int nInd = 0;
  for(int i=1;i<=n;i++){
    for(int j=1;j<=m;j++){
      nInd = (i-1)*m+(j-1);

      if(nModel == 1){
        eta[nInd]=_mu[nInd];
      }else{

        eta[nInd] = log(_mu[nInd]);
      }


    }
  }
  return;
}


void Dist::Setdeta_dmu(double* _mu){

  int nInd = 0;

  for(int i=1;i<=n;i++){
    for(int j=1;j<=m;j++){
      nInd = (i-1)*m+(j-1);

      if(nModel == 1){
        deta_dmu[nInd] = 1;
      }else{
        deta_dmu[nInd] = 1/_mu[nInd];
      }

    }
  }
  return;
}

void Dist::Sets(double* _mu, double* _y){

  ///// Normal Case
  //int nInd=0;
  int nInd = 0;

  if(nModel == 1){
    for(int i=1;i<=n;i++){
      for(int j=1;j<=m;j++){
        nInd = (i-1)*m+(j-1);

        if(nModel == 1){
          s[nInd] = _y[nInd];
        }else{
          s[nInd] = eta[nInd] + (_y[nInd] - _mu[nInd])*deta_dmu[nInd];
        }


      }
    }
  }


  return;
}



double* Dist::GetW(){
  return W;
}

double* Dist::Geteta(){
  return eta;
}

double* Dist::Getdeta_dmu(){
  return deta_dmu;
}

double* Dist::Gets(){
  return s;
}









