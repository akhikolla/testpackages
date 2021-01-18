
#include "Dist.h"
#include "Dist2.h"
#include "FuncLib.h"


mat GetWtilde(double* W, double* lambda, int n, int m){


  mat Wtilde(n,1);

  double tmp=0;
  int nInd=0;

  for(int i=1;i<=n;i++){
    tmp = 0;
    for(int j=1;j<=m; j++){
      nInd = (i-1)*m+(j-1);
      tmp += W[nInd];
    }

    Wtilde(i-1,0) = 1/(lambda[i-1])+tmp;

  }

  return Wtilde;
}


mat GetXWs(mat X, double* W, double* s){

  int nm = X.n_rows;
  int p = X.n_cols;
  int nInd = 0;


  //mat W(nm, nm, fill::zeros);
  mat Ws(nm, 1, fill::zeros);

  mat out(p, 1, fill::zeros);

  for(int i=1;i<=nm; i++){
    //W(i-1,i-1) = _W[i-1];
    Ws(i-1,0) = s[i-1] *  W[i-1];

  }

  out = X.t()*Ws;
  return out;
}

mat GetZWs(double* W, double* s, double* psi, double* lambda, int n, int m){

  mat out(n, 1, fill::zeros);

  int nInd = 0;
  double tmp = 0;

  for(int i=1;i<=n;i++){

    tmp = 0;
    for(int j=1;j<=m;j++){
      nInd = (i-1)*m+(j-1);
      tmp += W[nInd] * s[nInd];
    }
    out(i-1,0) = tmp + psi[i-1] / lambda[i-1];
  }

  return out;

}


mat GetZWX(mat X, double* _W, int n, int m){

  int p = X.n_cols;
  int nInd = 0;

  mat out(n, p, fill::zeros);
  mat Xi(m, p, fill::zeros);
  mat wi(1, m, fill::zeros);

  for(int i=1;i<=n;i++){

    Xi = X.submat((i-1)*m, 0, i*m-1 ,(p-1));

    for(int j=1;j<=m; j++){
      wi(0,j-1) = _W[(i-1)*m+(j-1)];
    }
    out.row(i-1) = wi*Xi;
  }

  return out;
}






bool CppGetDiag(double* out, double* _x, double* _y, double* _mu, double* lambda, int* _n, int* _m, int* _p){

  int n =  _n[0];
  int m =  _m[0];
  int p = _p[0];
  int np = n+p;


  //int nModel =_nModel[0];
  //int nModel2 = _nModel2[0];

  int nModel = 1;
  int nModel2 = 1;

  int nm = n*m;
  int nInd = 0;

  /////////////////////////
  double dummyW=0;
  double dummys=0;
  double dummydet=0;
  double dummyeta=0;
  double dummypsi=0;
  // /////////////////////////////

  double* W = &dummyW;

  double* eta = &dummys;
  double* deta_dmu = &dummydet ;
  double* s = &dummyeta;
  double* psi = &dummypsi;

  /////////////////////////////// Set up Main distribution and the second distribution

  Dist dist(n,m, nModel);

  dist.SetW(_mu);
  W = dist.GetW();

  dist.Seteta(_mu);
  eta = dist.Geteta();

  dist.Setdeta_dmu(_mu);
  deta_dmu = dist.Getdeta_dmu();

  dist.Sets(_mu, _y);
  s = dist.Gets();
  ////////////////////////////////////////

  Dist2 dist2(n, m, nModel2);

  dist2.Setpsi();
  psi = dist2.Getpsi();


  ///////////////////////////////////

  /////////// Get Wtilde

  mat Wtilde(n,1, fill::zeros);

  Wtilde = GetWtilde(W, lambda, n, m);

  double tmp=0;


  mat X(n*m,p, fill::zeros);

  int nInd2 = 0;

  for(int i=1;i<=n;i++){
    for(int j=1;j<=m;j++){
      nInd = (i-1)*m+(j-1);

      for(int k=1;k<=p;k++){
        nInd2 = (i-1)*m*p+(j-1)*p+(k-1);
        X( nInd, k-1) = _x[nInd2];
      }
    }
  }


  mat XWX(p,p, fill::zeros);

  /////////////////// Get XWX
  XWX = X.t()*X;
  //////////////////////

  mat XWs(p,1, fill::zeros);
  mat ZWs(n,1, fill::zeros);
  mat ZWX(n,p, fill::zeros);

  XWs = GetXWs(X, W, s);
  ZWs = GetZWs(W, s, psi, lambda, n, m);
  ZWX = GetZWX(X, W, n, m);

  mat onem(1,m, fill::zeros);

  mat onem2(m,1);
  onem2.ones();

  mat F11(p,p, fill::zeros);
  F11 = XWX;

  mat A12(p,n, fill::zeros );
  mat F12(p,n, fill::zeros);

  mat A21(n,p, fill::zeros);
  mat F21(n,p, fill::zeros);


  mat Xi(m,p, fill::zeros);

  mat F22(n,1, fill::zeros);

  //  double tmp=0;

  for(int i=1;i<=n;i++){
    //    tmp = 0;

    F22(i-1,0) = 1/Wtilde(i-1,0);
    Xi = X.submat( (i-1)*m, 0, i*m-1, (p-1) );

    for(int j=1; j<=m; j++){
      nInd = (i-1)*m+(j-1);
      onem(0,j-1) = W[nInd];
    }

    A21.row(i-1) = onem*Xi*F22(i-1,0);
    A12.col(i-1) = Xi.t()*onem.t();

    out[i-1] = F22(i-1,0);

  }

  double factor_ij = 0;

  for(int i=n;i>0;i--){

    for(int j=p;j>0;j--){

      factor_ij = A12(j-1,i-1);
      F12(j-1,i-1) = -factor_ij*F22(i-1,0);
      XWX.row(j-1) -= factor_ij*A21.row(i-1);
    }

  }
  //mat Idn(p,p);
  //Idn.eye();
  //F11 = solve(XWX, Idn); //  inv(XWX);
  F11 = inv(XWX);
  F12 = F11*F12;
  F21 = F12.t();

  mat bhat(p,1, fill::zeros);
  mat vhat(n,1, fill::zeros);
  mat tmpSecond(n,1, fill::zeros);

  bhat = F11*XWs + F12*ZWs;

  for(int i=1;i<=p;i++){
    out[i-1] = bhat(i-1,0);
  }

  tmpSecond = ZWs - ZWX*bhat;

  tmp=0;

  for(int i=1;i<=n;i++){
    tmp = 0;
    for(int j=1;j<=m;j++){
      tmp += W[(i-1)*m+(j-1)];
    }
    tmp += 1/lambda[i-1];
    vhat(i-1,0) = tmpSecond(i-1,0)/tmp;
    out[(i-1)+p] = vhat(i-1,0);
  }


  for(int j=p;j>0;j--){

    for(int i=n;i>0;i--){
      factor_ij = A21(i-1,j-1);
      F22(i-1,0) -= factor_ij* F12(j-1, i-1);

    }
  }



  //double* out;

  mat xijp(1,p, fill::zeros);
  mat f21r(1,p, fill::zeros);

  mat tmpmat(1,1);
  mat tmpmat2(1,1);

  //double tmp =0;
  tmp=0;
  double tmp2 =0;

  int nIndij = 0;

  for(int i=1; i<=n; i++){

    Xi.zeros();
    Xi = X.submat( (i-1)*m, 0, i*m-1, (p-1) );

    f21r.zeros();
    f21r = F21.row(i-1);

    for(int j=1; j<=m; j++){

      xijp.zeros();
      xijp = Xi.row(j-1);

      tmpmat2.zeros();
      tmpmat2 = 2*f21r * xijp.t();

      tmpmat.zeros();
      tmpmat = xijp * F11 * xijp.t();

      tmp=0;
      tmp = tmpmat(0,0);

      tmp2=0;
      tmp2 = tmpmat2(0,0);

      nIndij = (i-1)*m+(j-1)+np;

      out[nIndij] = ( tmp +  F22(i-1,0) + tmp2 ) ;
    }
    out[n*m+(i-1)+np] = F22(i-1,0) / lambda[i-1];
  }


  return 1;





}
