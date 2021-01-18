#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;
#include "penaliseTEMP.h"
Rcpp::NumericVector gFc(Rcpp::List liste, int l) {
  int lg = liste.size();
  if (l > lg) {
    Rcpp::stop("List index out of bounds");
  }
  l = l - 1;
  Rcpp::NumericVector v(liste[l]);
  return(v);
}//End gFc

double my_fct(double ro, void *params){
  struct my_fct_params *p = (struct my_fct_params *) params;
  int n = p->n;
  VectorXd Z1 = p->Z1;
  NumericVector d = p->d;
  NumericVector sqd = p->sqd;
  MatrixXd Q = p->Q;
  MatrixXd kv = p->kv;
  double muv = p->muv;
  
  NumericVector vp(n);
  vp = (muv*muv)/d + ro;
  
  NumericVector indvp;
  indvp = 1/(vp * sqd);
  
  NumericMatrix dg; dg = diag(indvp);
  
  VectorXd xChap;
  xChap = muv * Q * as<MatrixXd>(dg) * Q.transpose() * Z1;
  
  double NxChap; NxChap = xChap.norm()-1;
  return(NxChap);
}//End my_fct
double uniroot(int n,VectorXd Z1, NumericVector d,NumericVector sqd , MatrixXd Q,
               MatrixXd kv , double muv, double t0, double t1){
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r =0;
  double x_lo = t0, x_hi = t1;
  gsl_function F;
  
  struct my_fct_params params = {n,Z1, d,sqd , Q, kv, muv};
  F.function = &my_fct;
  F.params = &params;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
    
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  return r;
}//End uniroot
SEXP optV(int n,VectorXd Z1, NumericVector d,NumericVector sqd ,MatrixXd Q
            ,MatrixXd kv , double muv, double gamav){
  List l;
  bool a = true ;
  
  double nrm;
  nrm = sqrt(Z1.transpose()*kv*Z1);
  double mv;
  mv = nrm/muv;
  double x0; x0 = pow(mv,2);
  if(x0<=1){return List::create(Named("a",a),Named("crit",R_NilValue));}
  else {
    double ro = 1;
    double NxChap;
    do{
      ro = 10 * ro;
      NumericVector vp;
      vp = pow(muv , 2)/d + ro;
      NumericVector indvp;
      indvp = 1/(vp * sqd);
      NumericMatrix dg; dg = diag(indvp);
      VectorXd xChap;
      xChap = muv * Q * as<MatrixXd>(dg) * Q.transpose() * Z1;
      
      NxChap = xChap.norm();
    } while (NxChap>=1);
    if(ro==10){
      double t0 = 0;
      double t1 = ro;
      double root;
      root = uniroot(n,Z1,d,sqd ,Q,kv ,muv,t0,t1);
      NumericVector vp;
      vp = (muv*muv) + (root*d);
      NumericVector indvp;
      NumericMatrix dg;
      indvp= 1/(vp);
      dg = diag(indvp);
      VectorXd Z1Chap;
      Z1Chap = muv*muv * Q * as<MatrixXd>(dg) * Q.transpose() * Z1;
      
      VectorXd difZ;
      difZ = Z1-Z1Chap;
      double scr; scr = difZ.norm();
      double SCRChap; SCRChap = pow(scr,2);
      if(SCRChap<=pow(gamav,2)){
        l = List::create(Named("a",a),Named("crit",SCRChap-pow(gamav,2)));
      }else{
        l = List::create(Named("a",!a),Named("crit",SCRChap-pow(gamav,2)));
      }
    }
    if(ro>10){
      double t0 = ro/10;
      double t1 = ro;
      double root;
      root = uniroot(n,Z1,d,sqd ,Q,kv ,muv,t0,t1);
      NumericVector vp(n);
      vp = (muv*muv) + (root*d);
      NumericVector indvp;
      indvp= 1/(vp);
      NumericMatrix dg;
      dg = diag(indvp);
      VectorXd Z1Chap;
      Z1Chap = muv*muv * Q * as<MatrixXd>(dg) * Q.transpose() * Z1;
      VectorXd difZ;
      difZ = Z1-Z1Chap;
      double scr; scr = difZ.norm();
      double SCRChap; SCRChap = pow(scr,2);
      if(SCRChap<=pow(gamav,2)){
        l = List::create(Named("a",a),Named("crit",SCRChap-pow(gamav,2)));
      }else{
        l = List::create(Named("a",!a),Named("crit",SCRChap-pow(gamav,2)));
      }
    }
  }
  return l;
}//End optV

int rvfct_f(const gsl_vector * x, void *params, gsl_vector * f){
  NumericVector d = ((struct fparams *) params)->d;
  MatrixXd Q = ((struct fparams *) params)->Q;
  MatrixXd kv = ((struct fparams *) params)->kv;
  VectorXd R = ((struct fparams *) params)->R;
  double gamav = ((struct fparams *) params)->gamav;
  double muv = ((struct fparams *) params)->muv;
  
  const double x1 = gsl_vector_get (x, 0);
  const double x2 = gsl_vector_get (x, 1);
  NumericVector vp;
  vp = (1+x1)*d+x2;
  NumericVector invp;
  invp = 1/vp;
  NumericMatrix D;
  D = diag(invp);
  
  VectorXd teta;
  teta = Q * as<MatrixXd>(D) * Q.transpose() * R;//r small
  VectorXd nm1;
  nm1 = kv * teta;
  double Norm1;
  Norm1 = nm1.norm();
  double Norm2;
  Norm2 = sqrt(teta.transpose()* kv *teta);
  const double y0 = x1*Norm1-gamav/2;
  const double y1 = x2*Norm2-muv/2;
  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  
  return GSL_SUCCESS;
}//End rvfct_f
int rvfct_df(const gsl_vector * x, void *params, gsl_matrix * J){
  NumericVector d = ((struct fparams *) params)->d;
  MatrixXd Q = ((struct fparams *) params)->Q;
  MatrixXd kv = ((struct fparams *) params)->kv;
  VectorXd R = ((struct fparams *) params)->R;
  
  const double x1 = gsl_vector_get (x, 0);
  const double x2 = gsl_vector_get (x, 1);
  NumericVector vp; vp = (1+x1)*d+x2;
  NumericVector invp; invp = 1/vp;
  NumericMatrix D; D = diag(invp);
  VectorXd teta;
  teta = Q * as<MatrixXd>(D) * Q.transpose() * R;// r small
  
  VectorXd v1;
  v1 = kv * teta;
  double Norm1;
  Norm1 = v1.norm();
  double Norm2;
  Norm2 = sqrt(teta.transpose()*v1);
  
  NumericVector invp2; invp2 = (1/vp)*(1/vp);
  NumericVector dinvp2; dinvp2 = d*invp2;
  NumericMatrix D2; D2 = diag(dinvp2);
  
  VectorXd dalphadx1;
  dalphadx1 = -(Q * as<MatrixXd>(D2) * Q.transpose()  * R);
  NumericMatrix D3; D3 = diag(invp2);
  
  VectorXd dalphadx2;
  dalphadx2 = -(Q * as<MatrixXd>(D3) * Q.transpose()  * R);
  
  NumericVector d2; d2 = d*d;
  NumericMatrix DD; DD = diag(d2);
  MatrixXd cK2; cK2 = Q * as<MatrixXd>(DD) * Q.transpose();
  VectorXd v2;
  v2 = cK2 * teta;
  double dNorm1dx1;
  dNorm1dx1 = (1/Norm1)*dalphadx1.transpose()*v2;
  double dNorm2dx1;
  dNorm2dx1 = (1/Norm2)*dalphadx1.transpose()*v1;
  double dNorm1dx2;
  dNorm1dx2 = (1/Norm1)*dalphadx2.transpose()*v2;
  double dNorm2dx2;
  dNorm2dx2 = (1/Norm2)*dalphadx2.transpose()*v1;
  const double df00 = Norm1 + x1*dNorm1dx1;
  const double df01 = x1*dNorm1dx2;
  const double df10 = x2*dNorm2dx1;
  const double df11 = Norm2 + x2*dNorm2dx2;
  gsl_matrix_set (J, 0, 0, df00);
  gsl_matrix_set (J, 0, 1, df01);
  gsl_matrix_set (J, 1, 0, df10);
  gsl_matrix_set (J, 1, 1, df11);
  
  return GSL_SUCCESS;
}//End rvfct_df
int rvfct_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J){
  rvfct_f (x, params, f);
  rvfct_df (x, params, J);
  
  return GSL_SUCCESS;
}//End rvfct_fdf
SEXP nleqslvgnewton(NumericVector xstart, NumericVector d, MatrixXd Q, MatrixXd kv,
                    VectorXd R, double gamav, double muv){
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  
  int status;
  size_t iter = 0;
  const size_t n = 2;
  struct fparams p = {d,Q,kv,R,gamav,muv};
  gsl_multiroot_function_fdf f = {&rvfct_f,
                                  &rvfct_df,
                                  &rvfct_fdf,
                                  n, &p};
  
  double x_init[2] = {xstart(0), xstart(1)};
  gsl_vector *x = gsl_vector_alloc (n);
  
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  
  T = gsl_multiroot_fdfsolver_gnewton;
  s = gsl_multiroot_fdfsolver_alloc (T, n);
  gsl_multiroot_fdfsolver_set (s, &f, x);
  double r1;
  double r2;
  do
  {
    iter++;
    
    status = gsl_multiroot_fdfsolver_iterate (s);
    r1 = gsl_vector_get (s->x, 0);
    r2 = gsl_vector_get (s->x, 1);
    if (status)
      break;
    
    status = gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 500);
  
  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free (x);
  
  NumericVector v = NumericVector::create(r1,r2);
  return List::create(Named("status",status),Named("x",v));
}//End nleqslvgnewton
SEXP nleqslvhybrids(NumericVector xstart, NumericVector d, MatrixXd Q, MatrixXd kv,
                    VectorXd R, double gamav, double muv){
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  
  int status;
  size_t iter = 0;
  const size_t n = 2;
  struct fparams p = {d,Q,kv,R,gamav,muv};
  gsl_multiroot_function f = {&rvfct_f, n, &p};
  
  double x_init[2] = {xstart(0), xstart(1)};
  gsl_vector *x = gsl_vector_alloc (n);
  
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);
  double r1;
  double r2;
  do
  {
    iter++;
    
    status = gsl_multiroot_fsolver_iterate (s);
    r1 = gsl_vector_get (s->x, 0);
    r2 = gsl_vector_get (s->x, 1);
    if (status)
      break;
    
    status = gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 500);
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  
  NumericVector v = NumericVector::create(r1,r2);
  return List::create(Named("status",status),Named("x",v));
}//End nleqslvhybrids
SEXP nleqslvbroyden(NumericVector xstart, NumericVector d, MatrixXd Q, MatrixXd kv,
                    VectorXd R, double gamav, double muv){
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  
  int status;
  size_t iter = 0;
  const size_t n = 2;
  struct fparams p = {d,Q,kv,R,gamav,muv};
  gsl_multiroot_function f = {&rvfct_f, n, &p};
  
  double x_init[2] = {xstart(0), xstart(1)};
  gsl_vector *x = gsl_vector_alloc (n);
  
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  
  T = gsl_multiroot_fsolver_broyden;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);
  double r1;
  double r2;
  do
  {
    iter++;
    
    status = gsl_multiroot_fsolver_iterate (s);
    r1 = gsl_vector_get (s->x, 0);
    r2 = gsl_vector_get (s->x, 1);
    if (status)
      break;
    
    status = gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 500);
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  
  NumericVector v = NumericVector::create(r1,r2);
  return List::create(Named("status",status),Named("x",v));
}//End nleqslvbroyden
SEXP resiv(int n,NumericVector d, MatrixXd Q,MatrixXd kv, VectorXd R
             , double gamav, double muv, NumericVector tetav){
  bool cvge = true;
  NumericVector xme(2);
  VectorXd ttv(n);
  ttv = as<VectorXd>(tetav);
  double nrm;
  nrm = ttv.norm();
  double hyp; hyp = (nrm*nrm)/n;
  double ophyp; ophyp = (1/pow(10,16));
  if(hyp < ophyp){
    tetav = tetav + Rcpp::rnorm(n);
  }
  VectorXd tetv(n);
  tetv = as<VectorXd>(tetav);
  VectorXd mnorm1(n);
  mnorm1 = kv*tetv;
  double Norm1;
  Norm1 = mnorm1.norm();
  double  Norm2;
  Norm2 = sqrt(tetv.transpose()* kv *tetv);
  double x1;
  x1 = gamav/(2*Norm1);
  double x2;
  x2 = muv/(2*Norm2);
  NumericVector xstart(2); xstart[0]=x1; xstart[1]=x2;
  // try nleqslv with Generalized Newton's method
  List sln1;
  sln1 = nleqslvgnewton(xstart, d, Q, kv, R, gamav, muv);
  int status; status = sln1["status"];
  NumericVector x; x = sln1["x"];
  if(status==0){
    return List::create(Named("x",x),Named("cvge",cvge));
  }// If no cvge with gNewton's method
  if(status!=0){
    // try nleqslv with Broyden method
    List sln2;
    sln2 = nleqslvbroyden(xstart, d, Q, kv, R, gamav, muv);
    
    int status; status = sln2["status"];
    NumericVector x; x = sln2["x"];
    if(status==0){
      return List::create(Named("x",x),Named("cvge",cvge));
    }//If no cvge with gNewton's and Broyden method
    if(status!=0){
      // try nleqslv with Hybrid method
      List sln3;
      sln3 = nleqslvhybrids(xstart, d, Q, kv, R, gamav, muv);
      
      int status; status = sln3["status"];
      NumericVector x; x = sln3["x"];
      if(status==0){
        return List::create(Named("x",x),Named("cvge",cvge));
      }//If no cvge with gNewton's and Broyden and Hybrid method
      if(status!=0){
        return List::create(Named("cvge",!cvge));
      }
    }
  }
  return 0;
}//End resiv
// [[Rcpp::export]]
SEXP penMetaMod_cpp(NumericVector Y,List matZ,List k_v,StringVector namG, List resg,
                    NumericVector gamma, NumericVector mu,NumericVector gama_v,
                    NumericVector mu_v,int maxIter,
                    bool verbose, bool calcStwo){
  int n; n = Y.size();
  int vMax; vMax = namG.size();
  double eps; eps = 0.0001;
  
  if((gama_v(0)&&mu_v(0))==0){
    gama_v = rep(1,vMax);
    mu_v = rep(1,vMax);
  }
  int ls; ls =gamma.size();
  int lmu; lmu = mu.size();
  List Lme(lmu*ls);
  int sme=0;
  for(int lm=0;lm<lmu;lm++ ){
    double mu_i = sqrt(n)*mu[lm];
    List grp; grp = resg[lm];
    NumericVector supp_in; supp_in = grp["supp"];
    
    
    NumericVector teta_gr(vMax*n);
    teta_gr = gFc(grp,2);
    teta_gr.attr("dim") = Dimension(vMax, n);
    NumericMatrix teta_in(vMax,n); teta_in = as<NumericMatrix>(teta_gr);
    
    NumericVector fitv_gr(n*vMax);
    fitv_gr = gFc(grp,3);
    fitv_gr.attr("dim") = Dimension(n, vMax);
    NumericMatrix fitv_in(n,vMax); fitv_in = as<NumericMatrix>(fitv_gr);
    
    NumericVector NormH_in(vMax);
    NormH_in = grp["Norm.H"];
    
    for(int il=0;il<ls;il++ ){
      double gama; gama = sqrt(n)*gamma[il];
      if(gama!=0){
        List resK;
        
        NumericVector muvs; muvs = mu_i*mu_v;
        NumericVector gamavs; gamavs = gama*gama_v;
        
        double c;
        double SCR;
        double crit;
        NumericVector R(n);
        VectorXd r(n);
        MatrixXd kv(n,n);
        VectorXd Z1(n);
        double muv;
        double gamav;
        List a;
        List sln;
        
        double nrmRel;
        double rDiffCrit;
        double nrmRel2;
        double rDiffCrit2;
        NumericMatrix old_teta(vMax,n);
        
        bool convergence;
        
        NumericVector zerosv(vMax);
        NumericVector zeros(n);
        NumericVector tetav(n);
        
        NumericMatrix teta(vMax,n);
        for(int l=0; l<vMax;l++){
          teta(l,_) = teta_in(l,_);
        }
        NumericMatrix fitv(n,vMax);
        for(int l=0; l<vMax;l++){
          fitv(_,l) = fitv_in(_,l);
        }
        NumericVector Norm_H(vMax);
        for(int l=0; l<vMax;l++){
          Norm_H(l) = NormH_in(l);
        }
        NumericVector Norm_n(vMax);
        Norm_n.names() = namG;
        Norm_H.names() = namG;
        
        //Step 1
        NumericVector supp;
        for(int s=0;s<supp_in.size();s++){
          supp.push_back(supp_in[s]);
        }
        double penalite = 0;
        NumericVector fitted(n);
        if(supp.size()>0){
          for(NumericVector::iterator v = supp.begin(); v != supp.end(); ++v) {
            resK = matZ[(*v-1)];
            NumericVector d;
            d = gFc(resK,1);
            MatrixXd Q; Q = resK["Q"];
            tetav = teta((*v-1),_);
            
            Norm_n(*v-1) = sqrt(sum(pow(fitv(_,(*v-1)),2)));
            fitted = fitted + fitv(_,(*v-1));
            penalite = penalite + gamavs(*v-1) * Norm_n(*v-1) + muvs(*v-1) * Norm_H(*v-1);
          }
        }
        c = mean(Y-fitted);
        NumericVector ycf; ycf = Y-c-fitted;
        VectorXd ycfX; ycfX = as<VectorXd>(ycf);
        double ycfnorm; ycfnorm = ycfX.norm();
        SCR = pow(ycfnorm,2);
        crit = SCR + penalite;
        int ii; ii=0;
        NumericVector support;
        do{
          convergence = true;
          ii = ii+1;
          if(verbose){Rcout<<"Iteration number in Step1 is "<<ii<<"\n";}
          for(int l=0; l<vMax;l++){
            old_teta(l,_) = teta(l,_);
          }
          double critii; critii = crit;
          if(supp.size()>0){
            for(NumericVector::iterator v = supp.begin(); v != supp.end(); ++v){
              resK = matZ[(*v-1)];
              NumericVector d(n);
              d = gFc(resK,1);
              MatrixXd Q; Q = resK["Q"];
              MatrixXd kv = k_v[(*v-1)];
              tetav = old_teta(*v-1,_);
              NumericVector S(n);
              int ss=0;
              if(supp.size()==2){
                if(supp[0] == *v){
                  ss = supp[1]-1;
                  S = fitv(_,ss);
                }
                if(supp[1] == *v){
                  ss = supp[0]-1;
                  S = fitv(_,ss);
                }
              }
              if(supp.size()>2){
                NumericVector sv;
                for(int is=0;is<supp.size();is++){
                  if(supp[is] != *v){
                    ss = supp[is];
                    sv.push_back(ss);
                  }
                }
                for(NumericVector::iterator ivv = sv.begin(); ivv != sv.end(); ++ivv){
                  S += fitv(_,(*ivv-1));
                }
              }
              if(supp.size()==1){
                S = zeros;
              }
              c = mean(Y-S-fitv(_,(*v-1)));
              R = Y-c-S;
              r = as<VectorXd>(R);
              Z1 = 2*r;
              NumericVector sqd = sqrt(d);
              muv = muvs[(*v-1)];
              gamav = gamavs[(*v-1)];
              
              a = optV(n,Z1,d,sqd ,Q,kv , muv,gamav);
              bool zero; zero = a["a"];
              if(zero){
                teta((*v-1),_) = zeros;
                fitv(_,(*v-1)) = zeros;
              }
              if (!zero) {
                sln = resiv(n,d, Q,kv, r, gamav, muv, tetav);
                bool cvge; cvge = sln["cvge"];
                if (!cvge) {
                  if(verbose){
                    Rcout<<"ii : "<<ii<<"convergence failed for group : "<<namG[(*v-1)]<<"\n";
                  }
                  teta((*v-1),_) = tetav;
                  convergence = false;
                }
                if (cvge) {
                  
                  NumericVector x; x = sln["x"];
                  NumericVector vp; vp = (1+x(0)) * d+x(1);
                  NumericVector invp; invp = 1/vp;
                  NumericMatrix ivp; ivp = diag(invp);
                  VectorXd vcc; vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
                  NumericVector vcc0; vcc0 = vcc;
                  teta((*v-1),_) = vcc0;
                  VectorXd ckvteta; ckvteta = kv * vcc;
                  NumericVector ckvteta0; ckvteta0 = ckvteta;
                  fitv(_,(*v-1)) = ckvteta0;
                }
              }// end if (!zero) {
            }// end for(v in supp){
          }// end if (length(supp) >0)
          //Rcout<<"convergence is"<<convergence<<"\n";
          Norm_H = 0*Norm_H;
          Norm_n = 0*Norm_n;
          if(convergence){
            MatrixXd ttdef; ttdef = as<MatrixXd>(teta)-as<MatrixXd>(old_teta);
            double nrm; nrm = ttdef.norm();
            double nrmii; nrmii = (as<MatrixXd>(old_teta)).norm();
            if (nrmii>0) {nrmRel = (nrm/nrmii)*(nrm/nrmii);}
            if (nrmii==0) {nrmRel = nrm*nrm;}
            penalite = 0;
            fitted = rep(c,n);
            NumericVector suppStep1;
            if(supp.size()>0){
              for(NumericVector::iterator v = supp.begin(); v != supp.end(); ++v){
                if(sum(pow(fitv(_,(*v-1)),2))>1/pow(10,8)){
                  suppStep1.push_back(*v);
                }
              }
            }
            support = suppStep1;
            if(suppStep1.size()>0){
              for(NumericVector::iterator v = suppStep1.begin(); v != suppStep1.end(); ++v){
                Norm_n(*v-1)=sqrt(sum(pow(fitv(_,(*v-1)),2)));
                kv = k_v[(*v-1)];
                NumericVector tetav; tetav = teta((*v-1),_);
                VectorXd tetv; tetv = as<VectorXd>(tetav);
                Norm_H(*v-1) = sqrt(tetv.transpose()*kv*tetv);
                fitted = fitted + fitv(_,(*v-1));
                penalite = penalite+gamavs[(*v-1)]*Norm_n(*v-1)+muvs[(*v-1)]*Norm_H(*v-1);
              }
            }
            NumericVector yf; yf = Y-fitted;
            int yfsize; yfsize = yf.size();
            yf.attr("dim") = Dimension(yfsize,1);
            MatrixXd yfm; yfm = as<MatrixXd>(yf);
            SCR = pow(yfm.norm(),2);
            crit = SCR+penalite;
            rDiffCrit = abs((critii-crit)/critii);
            if(verbose){
              Rcout<<"Active support in Step1 is "<<support<<"\n";
              Rcout<<"rDiffCrit=crit_i-crit_(i-1)/crit_i in Step1 is "<<rDiffCrit<<"\n";
              Rcout<<"nrmRel=(norm(teta_i-teta_(i-1))/norm(teta_(i-1)))^2 in Step1 is "<<nrmRel<<"\n";
            }
            // 0 is false and 1 is true : 0&1 = 0 ,...
            int brk; brk = ((nrmRel<eps)&(rDiffCrit<eps));
            // 0&&1 = false
            if (brk&&(ii>=1)) break;
          }// end  if (convergence){
        } while (ii<=maxIter);
        if (ii>maxIter) {
          convergence = false;
        }
        //  End of Step 1
        
        int iiStep1 = ii;
        double nrmRelStep1 = nrmRel;
        double rDiffCritStep1 = rDiffCrit;
        if(!calcStwo){
          if(convergence){
            colnames(fitv) = namG;
            StringVector Nsupp;
            if(support.size()>0) Nsupp = namG[(support-1)];
            rownames(teta) = namG;
            List iter = List::create(Named("maxIter",maxIter),Named("Step1",iiStep1));
            List rtrn = List::create(Named("intercept",c),Named("teta",teta),Named("fit.v",fitv),
                                           Named("fitted",fitted),Named("Norm.n",Norm_n),Named("Norm.H",Norm_H),
                                           Named("supp",support),Named("Nsupp",Nsupp),Named("SCR",SCR),
                                           Named("crit",crit),Named("gamma.v",gamavs),Named("mu.v",muvs),
                                           Named("iter",iter),Named("convergence",convergence),
                                           Named("RelDiffCrit",rDiffCritStep1),Named("RelDiffPar",nrmRelStep1));
            Lme[il+lm+sme]=List::create(Named("mu",mu[lm]/sqrt(n)),Named("gamma",gamma[il]),Named("Meta-Model",rtrn));
          }//if(convergence)
          if(!convergence){
            StringVector Nsupp;
            if(support.size()>0) Nsupp = namG[(support-1)];
            rownames(teta) = namG;
            List iter = List::create(Named("maxIter",maxIter),Named("Step1",iiStep1));
            List rtrn = List::create(Named("intercept",c),Named("teta",teta),Named("supp",support),
                                           Named("Nsupp",Nsupp),Named("gamma.v",gamavs),Named("mu.v",muvs),
                                           Named("iter",iter),Named("convergence",convergence));
            Lme[il+lm+sme]=List::create(Named("mu",mu[lm]/sqrt(n)),Named("gamma",gamma[il]),Named("Meta-Model",rtrn));
          }//if(!convergence)
        }//if(!calcStwo)
        //  Step 2 :
        if(calcStwo){
          if(verbose){Rcout<<"Start step2 and supp of Step1 is "<<support<<"\n";}
          supp = zerosv;
          NumericVector supp2;
          ii = 0;
          double crit2; crit2 = crit;
          do{
            convergence = true;
            ii = ii+1;
            if(verbose){Rcout<<"Iteration number in Step2 is "<<ii<<"\n";}
            for(int l=0; l<vMax;l++){
              old_teta(l,_) = teta(l,_);
            }
            double critii; critii = crit2;
            for(int v=0; v<vMax; ++v){
              resK = matZ[v];
              NumericVector d;
              d = gFc(resK,1);
              NumericVector sqd(n);
              sqd = sqrt(d);
              MatrixXd Q; Q = resK["Q"];
              tetav = old_teta(v,_);
              MatrixXd kv = k_v[v];
              NumericVector sv(0);
              NumericVector S(n);
              for(int is=0;is<vMax;is++){
                if(is != v){
                  sv.push_back(is);
                }
              }
              for(NumericVector::iterator ivv = sv.begin(); ivv != sv.end(); ++ivv){
                S += fitv(_,*ivv);
              }
              c = mean(Y-S-fitv(_,v));
              R = Y-c-S;
              r = as<VectorXd>(R);
              Z1 = 2*r;
              muv = muvs[v];
              gamav = gamavs[v];
              a = optV(n,Z1,d,sqd ,Q,kv , muv,gamav);
              bool zero; zero =a["a"];
              
              if(zero){
                teta(v,_) = zeros;
                fitv(_,v) = zeros;
                supp[v]=0;
              }
              
              if (!zero) {
                supp[v]=1;
                sln = resiv(n,d, Q,kv, r, gamav, muv, tetav);
                bool cvge; cvge = sln["cvge"];
                if (!cvge) {
                  if(verbose){
                    Rcout<<"ii : "<<ii<<"convergence failed for group : "<<namG[v]<<"\n";
                  }
                  teta(v,_) = tetav;
                  convergence = false;
                }
                if (cvge) {
                  NumericVector x; x = sln["x"];
                  NumericVector vp; vp = (1+x(0)) * d+x(1);
                  NumericVector invp; invp = 1/vp;
                  NumericMatrix ivp; ivp = diag(invp);
                  VectorXd vcc; vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
                  NumericVector vcc0; vcc0 = vcc;
                  teta(v,_) = vcc0;
                  VectorXd ckvteta;
                  ckvteta = kv * vcc;
                  NumericVector ckvteta0; ckvteta0 = ckvteta;
                  fitv(_,v) = ckvteta0;
                }
              }// end if (!zero) {
            }//end for(v in 1:vMax) {
            Norm_H = 0*Norm_H;
            Norm_n = 0*Norm_n;
            if(convergence){
              NumericVector sq; sq = seq(1,vMax);
              supp2 = sq[(supp==1)];
              MatrixXd adef; adef = as<MatrixXd>(teta)-as<MatrixXd>(old_teta);
              double nrm; nrm = adef.norm();
              double nrmii; nrmii = (as<MatrixXd>(old_teta)).norm();
              if (nrmii>0) {nrmRel2 = (nrm/nrmii)*(nrm/nrmii);}
              if (nrmii==0) {nrmRel2 = nrm*nrm;}
              penalite = 0;
              fitted = rep(c,n);
              if(supp2.size()>0){
                for(NumericVector::iterator v = supp2.begin(); v != supp2.end(); ++v){
                  resK = matZ[(*v-1)];
                  NumericVector d;
                  d = gFc(resK,1);
                  MatrixXd Q; Q = resK["Q"];
                  MatrixXd kv = k_v[(*v-1)];
                  Norm_n(*v-1) = sqrt(sum(pow(fitv(_,(*v-1)),2)));
                  NumericVector tetav = teta((*v-1),_);
                  VectorXd tetv = as<VectorXd>(tetav);
                  double norm_hv = sqrt(tetv.transpose()*kv*tetv);
                  Norm_H(*v-1) = norm_hv;
                  fitted = fitted + fitv(_,(*v-1));
                  penalite = penalite+gamavs[(*v-1)]*Norm_n(*v-1)+muvs[(*v-1)]*Norm_H(*v-1);
                }
              }
              NumericVector yf; yf = Y-fitted;
              int yfsize; yfsize = yf.size();
              yf.attr("dim") = Dimension(yfsize,1);
              MatrixXd yfm; yfm = as<MatrixXd>(yf);
              SCR = pow(yfm.norm(),2);
              crit2 = SCR+penalite;
              rDiffCrit2 = abs((critii-crit2)/critii);
              if(verbose){
                Rcout<<"Active support in Step2 is "<<supp2<<"\n";
                Rcout<<"rDiffCrit=crit_i-crit_(i-1)/crit_i in Step2 is "<<rDiffCrit2<<"\n";
                Rcout<<"nrmRel=(norm(teta_i-teta_(i-1))/norm(teta_(i-1)))^2 in Step2 is "<<nrmRel2<<"\n";
              }
              // 0 is false and 1 is true : 0&1 = 0 ,...
              int brk; brk = ((nrmRel2<eps)&(rDiffCrit2<eps));
              // 0&&1 = false
              if (brk&&(ii>=1)) break;
            }// end  if (convergence) {
          } while (ii<=maxIter);
          if (ii>maxIter) {
            convergence = false;
          }
          //end of do
          double rDiffCritStep2 = rDiffCrit2;
          double nrmRelStep2 = nrmRel2;
          //  End of Step 2
          if(convergence){
            colnames(fitv) = namG;
            StringVector Nsupp;
            if(supp.size()>0) Nsupp = namG[(supp2-1)];
            rownames(teta) = namG;
            
            List RelDiffCrit = List::create(Named("Step1",rDiffCritStep1),Named("Step2",rDiffCritStep2));
            List RelDiffPar = List::create(Named("Step1",nrmRelStep1),Named("Step2",nrmRelStep2));
            List iter = List::create(Named("maxIter",maxIter),Named("Step1",iiStep1),Named("Step2",ii));
            List rtrn = List::create(Named("intercept",c),Named("teta",teta),Named("fit.v",fitv),
                                           Named("fitted",fitted),Named("Norm.n",Norm_n),Named("Norm.H",Norm_H),
                                           Named("supp",supp2),Named("Nsupp",Nsupp),Named("SCR",SCR),
                                           Named("crit",crit),Named("gamma.v",gamavs),Named("mu.v",muvs),
                                           Named("iter",iter),Named("convergence",convergence),
                                           Named("RelDiffCrit",RelDiffCrit),Named("RelDiffPar",RelDiffPar));
            Lme[il+lm+sme]=List::create(Named("mu",mu[lm]/sqrt(n)),Named("gamma",gamma[il]),Named("Meta-Model",rtrn));
          }//if(convergence)
          if(!convergence){
            StringVector Nsupp;
            if(supp.size()>0) Nsupp = namG[(supp2-1)];
            rownames(teta) = namG;
            List iter = List::create(Named("maxIter",maxIter),Named("Step1",iiStep1),Named("Step2",ii));
            List rtrn = List::create(Named("intercept",c),Named("teta",teta),Named("supp",supp2),
                                           Named("Nsupp",Nsupp),Named("gamma.v",gamavs),Named("mu.v",muvs),
                                           Named("iter",iter),Named("convergence",convergence));
            Lme[il+lm+sme]=List::create(Named("mu",mu[lm]/sqrt(n)),Named("gamma",gamma[il]),Named("Meta-Model",rtrn));
          }//if(!convergence)
        }//if(calcStwo)
      }else{
        Lme[il+lm+sme]=List::create(Named("mu",mu[lm]),Named("gamma",gamma[il]),Named("Meta-Model",resg[lm]));
      }
    }//End of for(all gamma)
    sme+=(ls-1);
  }//End of for(all mu)
  return Lme;
}//End penMetaMod_cpp
