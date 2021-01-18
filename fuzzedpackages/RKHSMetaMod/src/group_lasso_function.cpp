#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;
using namespace std;
struct tetaf_params{
  NumericVector d;
  NumericVector sqd;
  MatrixXd Q;
  VectorXd R;
  double mqn;
};
double tetaf (double ro, void *params){
  struct tetaf_params *p = (struct tetaf_params *) params;

  NumericVector d = p->d;
  NumericVector sqd = p->sqd;
  MatrixXd Q = p->Q;
  VectorXd R = p->R;
  double mqn = p->mqn;


  NumericVector vp; vp = (sqd*ro)/(d+ro);
  NumericMatrix D; D = diag(vp);
  VectorXd sqkvtetav_ro; sqkvtetav_ro = Q*as<MatrixXd>(D)*Q.transpose()*R;

  double y; y = (2*sqkvtetav_ro.norm())-mqn;
  return y;
}//End tetaf
SEXP slv(NumericVector d,NumericVector sqd, MatrixXd Q, VectorXd R, double mqn,
         double t0, double t1){
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r =0;
  double x_lo = t0, x_hi = t1;
  gsl_function F;

  struct tetaf_params params = {d,sqd, Q, R, mqn};
  F.function = &tetaf;
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
  return List::create(Named("status",status),Named("x",r));
}//End slv
// [[Rcpp::export]]
SEXP RKHSgrplasso(NumericVector Y, List Kv, double mu, int maxIter=1000, bool verbose=false){
  double eps=1e-4;
  int n = Y.size();
  List matZ; matZ = Kv[0];
  List resK(2);
  List L(10);
  StringVector namesGrp; namesGrp = Kv[1];
  int vMax = namesGrp.size();
  NumericMatrix teta(vMax, n);
  NumericMatrix fitv(n,vMax);
  NumericVector zeros(n);

  List k_v(vMax);
  NumericVector sqd(n);
  double crit;
  NumericVector yf; yf = Y-mean(Y);
  crit = pow(as<VectorXd>(yf).norm(),2)/sqrt(n);

  double nrmRel;
  double rDiffCrit;

  bool convergence;
  double c;
  NumericVector R(n);//Residuals
  VectorXd r(n);
  NumericMatrix old_teta(teta.rows(),teta.cols());
  //Define parameters
  NumericVector d(n);
  MatrixXd Q(n,n);
  MatrixXd kv(n,n);
  NumericVector tetav(n);
  VectorXd kvtetav(n);
  NumericVector kvtetav0(n);
  double mqn = mu*sqrt(n);

  double cnd;
  NumericVector abstv(n);
  double somm;
  NumericVector vp_init(n);
  NumericMatrix D_init(n,n);
  VectorXd sqkvtetav_ro_init(n);
  double tetaf_init;
  double ro;
  double tetaf;
  NumericVector vp(n);
  NumericMatrix D(n,n);
  VectorXd sqkvtetav_ro(n);
  double t0;
  double t1;
  List calcslv(2);
  int status;
  NumericMatrix ivp(n,n);
  VectorXd vcc(n);
  NumericVector vcc0(n);
  VectorXd kvteta(n);
  NumericVector kvteta0(n);
  double root;
  VectorXd sqkvtetav(n);
  double Norm;
  double ro_init;
  double norm_hv;
  NumericVector Norm_H(vMax);
  
  NumericVector fitted(n);
  double SCR;
  NumericVector support;
  if(mu==0){
    Rcpp::stop("mu should be positive.");
  }
  for(int v=0; v<vMax;v++){
    resK = matZ[v];
    d = resK["Evalues"];
    Q = resK["Q"];
    sqd = sqrt(d);
    NumericMatrix D; D = diag(d);
    kv = Q * as<MatrixXd>(D) * Q.transpose();
    k_v[v] = kv;
    tetav = old_teta(v,_);

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
    c = mean(Y-S-fitv(_,v));//estimation of f_0
    R = Y-c-S;//Residuals
    r=as<VectorXd>(R);
    cnd = 2*sqrt(r.transpose()*kv*r);
    if(cnd<=mqn){
      teta(v,_) = zeros;
      fitv(_,v) = zeros;
    }
    else{
      abstv = abs(tetav);
      somm = sum(abstv);
      if(somm==0){
        vp_init = (sqd)/(d+1);
        D_init = diag(vp_init);
        sqkvtetav_ro_init = Q*as<MatrixXd>(D_init)*Q.transpose()*r;
        tetaf_init = 2*sqkvtetav_ro_init.norm();
        if(tetaf_init<mqn){
          ro = 1;
          do{
            ro = 10 * ro;
            vp = (sqd*ro)/(d+ro);
            D = diag(vp);
            sqkvtetav_ro = Q*as<MatrixXd>(D)*Q.transpose()*r;
            tetaf = 2*sqkvtetav_ro.norm();
          } while (tetaf<mqn);
          t0 = ro/10;
          t1 = ro;
          calcslv = slv(d,sqd, Q, r, mqn, t0, t1);
          status = calcslv["status"];

          if(status==0){
            root = calcslv["x"];
            vp = 1/(d+root);
            ivp = diag(vp);
            vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
            vcc0 = vcc;
            teta(v,_) = vcc0;
            kvteta = kv * vcc;
            kvteta0 = kvteta;
            fitv(_,v) = kvteta0;
            }else{
            teta(v,_) = zeros;
            fitv(_,v) = zeros;
          }
        }//end if tetaf_init<mu*sqrt(n)

        if(tetaf_init>=mqn){
          t0 = 0;
          t1 = 1;
          calcslv = slv(d,sqd, Q, r, mqn, t0, t1);
          status = calcslv["status"];
          if(status==0){
            root = calcslv["x"];
            vp = 1/(d+root);
            ivp = diag(vp);
            vcc = Q * as<MatrixXd>(ivp) * Q.transpose() *r;
            vcc0 = vcc;
            teta(v,_) = vcc0;
            kvteta = kv * vcc;
            kvteta0 = kvteta;
            fitv(_,v) = kvteta0;
            }else{
            teta(v,_) = zeros;
            fitv(_,v) = zeros;
          }
        }//end if tetaf_init>=mu/2
      }//end if(somm==0)
      if(somm!=0){
        VectorXd tetv;
        tetv = as<VectorXd>(tetav);
        Norm = sqrt(tetv.transpose()* kv*tetv);
        ro_init = (mqn)/(2*(Norm));
        vp_init = (sqd*ro_init)/(d+ro_init);
        D_init = diag(vp_init);
        sqkvtetav_ro_init = Q*as<MatrixXd>(D_init)*Q.transpose()*r;

        tetaf_init = 2*sqkvtetav_ro_init.norm();
        if(tetaf_init<mqn){
          ro = ro_init;
          do{
            ro = 10 * ro;
            vp = (sqd*ro)/(d+ro);
            D = diag(vp);
            sqkvtetav_ro = Q*as<MatrixXd>(D)*Q.transpose()*r;

            tetaf = 2*sqkvtetav_ro.norm();
          } while (tetaf<mqn);
          t0 = ro/10;
          t1 = ro;
          calcslv = slv(d,sqd, Q, r, mqn, t0, t1);
          status = calcslv["status"];

          if(status==0){
            root = calcslv["x"];
            vp = 1/(d+root);
            ivp = diag(vp);
            vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
            vcc0 = vcc;
            teta(v,_) = vcc0;
            kvteta = kv * vcc;
            kvteta0 = kvteta;
            fitv(_,v) = kvteta0;
            }else{
            teta(v,_) = old_teta(v,_);
            kvtetav = kv*as<VectorXd>(tetav);
            kvtetav0 = kvtetav;
            fitv(_,v) = kvtetav0;
          }
        }//end if tetaf_init<mu*sqrt(n)

        if(tetaf_init>=mqn){
          ro = ro_init;
          do{
            ro = ro/10;
            vp = (sqd*ro)/(d+ro);
            D = diag(vp);
            sqkvtetav_ro = Q*as<MatrixXd>(D)*Q.transpose()*r;

            tetaf = 2*sqkvtetav_ro.norm();
          } while (tetaf>=mqn);
          t0 = ro;
          t1 = 10*ro;
          calcslv = slv(d,sqd, Q, r, mqn, t0, t1);
          status = calcslv["status"];

          if(status==0){
            root = calcslv["x"];
            vp = 1/(d+root);
            ivp = diag(vp);
            vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
            vcc0 = vcc;
            teta(v,_) = vcc0;
            kvteta = kv * vcc;
            kvteta0 = kvteta;
            fitv(_,v) = kvteta0;
            }else{
            teta(v,_) = old_teta(v,_);
            kvtetav = kv*as<VectorXd>(tetav);
            kvtetav0 = kvtetav;
            fitv(_,v) = kvtetav0;
          }
        }//end if tetaf_init>=mu*sqrt(n)
      }//end if(somm!=0)
    }//End else
  }//End v in 1:vMax
  //Second STEP
  NumericVector active;
  for(int vm=1;vm<vMax+1;vm++){
    NumericVector abstt(n); abstt = abs(teta((vm-1),_));
    double SuM; SuM = sum(abstt);
    if(SuM!=0){
      active.push_back(vm);
    }
  }
  if(active.size()==0){
    Rcout<<"All groups are null"<<"\n";
    colnames(fitv) = namesGrp;
    rownames(teta) = namesGrp;
    Norm_H.names() = namesGrp;
    L = List::create(Named("intercept",c),Named("teta",teta),Named("fit.v",fitv)
                       ,Named("crit",crit),Named("Norm.H",Norm_H)
                       ,Named("convergence",true),Named("MaxIter",1)
                       ,Named("nrmRel",0), Named("rDiffCrit",0));

  }else{
    int i=0;
    do{
      convergence = true;
      i = i+1;
      if(verbose){Rcout<<"Iteration number is "<<i<<"\n";}
      for(int l=0; l<vMax;l++){
        old_teta(l,_) = teta(l,_);
      }
      double criti; criti = crit;
      for(int v=0; v<vMax;v++){
        resK = matZ[v];
        d = resK["Evalues"];
        Q = resK["Q"];
        NumericVector sqdd = sqrt(d);
        tetav = old_teta(v,_);
        kv = k_v[v];
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
        c = mean(Y-S-fitv(_,v));//estimation of f_0
        R = Y-c-S;//Residuals
        r=as<VectorXd>(R);
        cnd = 2*sqrt(r.transpose()*kv*r);
        if(cnd<=mqn){
          teta(v,_) = zeros;
          fitv(_,v) = zeros;
        }
        else{
          abstv = abs(tetav);
          somm = sum(abstv);
          if(somm==0){
            vp_init = sqdd/(d+1);
            D_init = diag(vp_init);
            sqkvtetav_ro_init = Q*as<MatrixXd>(D_init)*Q.transpose()*r;

            tetaf_init = 2*sqkvtetav_ro_init.norm();
            if(tetaf_init<mqn){
              ro = 1;
              do{
                ro = 10 * ro;
                vp = (sqdd*ro)/(d+ro);
                D = diag(vp);
                sqkvtetav_ro = Q*as<MatrixXd>(D)*Q.transpose()*r;

                tetaf = 2*sqkvtetav_ro.norm();
              } while (tetaf<mqn);
              t0 = ro/10;
              t1 = ro;
              calcslv = slv(d,sqdd, Q, r, mqn, t0, t1);
              status = calcslv["status"];

              if(status==0){
                root = calcslv["x"];
                vp = 1/(d+root);
                ivp = diag(vp);
                vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
                vcc0 = vcc;
                teta(v,_) = vcc0;
                kvteta = kv * vcc;
                kvteta0 = kvteta;
                fitv(_,v) = kvteta0;
                }else{
                teta(v,_) = zeros;
                fitv(_,v) = zeros;
                convergence = false;
              }
            }//end if tetaf_init<mu*sqrt(n)

            if(tetaf_init>=mqn){

              t0 = 0;
              t1 = 1;
              calcslv = slv(d,sqdd, Q, r, mqn, t0, t1);
              status = calcslv["status"];

              if(status==0){
                root = calcslv["x"];
                vp = 1/(d+root);
                ivp = diag(vp);
                vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
                vcc0 = vcc;
                teta(v,_) = vcc0;
                kvteta = kv * vcc;
                kvteta0 = kvteta;
                fitv(_,v) = kvteta0;
                }else{
                teta(v,_) = zeros;
                fitv(_,v) = zeros;
                convergence = false;
              }
            }//end if tetaf_init>=mu*sqrt(n)
          }//end if(somm==0)
          if(somm!=0){
            VectorXd tetv;
            tetv=as<VectorXd>(tetav);
            Norm = sqrt(tetv.transpose()* kv *tetv);
            ro_init = (mqn)/(2*(Norm));
            vp_init = (sqdd * ro_init)/(d+ro_init);
            D_init = diag(vp_init);
            sqkvtetav_ro_init = Q*as<MatrixXd>(D_init)*Q.transpose()*r;

            tetaf_init = 2*sqkvtetav_ro_init.norm();
            if(tetaf_init<mqn){

              ro = ro_init;
              do{
                ro = 10 * ro;
                vp = (sqdd*ro)/(d+ro);
                D = diag(vp);
                sqkvtetav_ro = Q*as<MatrixXd>(D)*Q.transpose()*r;

                tetaf = 2*sqkvtetav_ro.norm();
              } while (tetaf<mqn);
              t0 = ro/10;
              t1 = ro;
              calcslv = slv(d,sqdd, Q, r, mqn, t0, t1);
              status = calcslv["status"];

              if(status==0){
                root = calcslv["x"];
                vp = 1/(d+root);
                ivp = diag(vp);
                vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
                vcc0 = vcc;
                teta(v,_) = vcc0;
                kvteta = kv * vcc;
                kvteta0 = kvteta;
                fitv(_,v) = kvteta0;
                }else{
                teta(v,_) = old_teta(v,_);
                kvtetav = kv*as<VectorXd>(tetav);
                kvtetav0 = kvtetav;
                fitv(_,v) = kvtetav0;
                convergence = false;
              }
            }//end if tetaf_init<mu*sqrt(n)

            if(tetaf_init>=mqn){
              ro = ro_init;
              do{
                ro = ro/10;
                vp = (sqdd*ro)/(d+ro);
                D = diag(vp);
                sqkvtetav_ro = Q*as<MatrixXd>(D)*Q.transpose()*r;

                tetaf = 2*sqkvtetav_ro.norm();

              } while (tetaf>=mqn);
              t0 = ro;
              t1 = 10*ro;
              calcslv = slv(d,sqdd, Q, r, mqn, t0, t1);
              status = calcslv["status"];

              if(status==0){
                root = calcslv["x"];
                vp = 1/(d+root);
                ivp = diag(vp);
                vcc = Q * as<MatrixXd>(ivp) * Q.transpose() * r;
                vcc0 = vcc;
                teta(v,_) = vcc0;
                kvteta = kv * vcc;
                kvteta0 = kvteta;
                fitv(_,v) = kvteta0;
                }else{
                teta(v,_) = old_teta(v,_);
                kvtetav = kv*as<VectorXd>(tetav);
                kvtetav0 = kvtetav;
                fitv(_,v) = kvtetav0;
                convergence = false;
              }
            }//end if tetaf_init>=mu*sqrt(n)
          }//end if(somm!=0)
        }//End else
      }//End v in vMax
      if(convergence){
        MatrixXd ttdif; ttdif = as<MatrixXd>(teta)-as<MatrixXd>(old_teta);
        double nrm; nrm = ttdif.norm();
        double nrmi; nrmi = (as<MatrixXd>(old_teta)).norm();
        if (nrmi==0) {
          nrmRel = nrm*nrm;
        }else{
          nrmRel = (nrm/nrmi)*(nrm/nrmi);
        }
        NumericVector supp;
        for(int v=0;v<vMax;v++){
          NumericVector abstt; abstt = abs(teta(v,_));
          double SuM; SuM = sum(abstt);
          if(SuM!=0){
            supp.push_back(v+1);
          }
        }
        double penalite = 0;
        fitted = rep(c,n);
        support = supp;
        NumericVector Norm_h(vMax);
        if(supp.size()>0){
          for(NumericVector::iterator v = supp.begin(); v != supp.end(); ++v){
            NumericVector tetav; tetav = teta((*v-1),_);
            VectorXd tetv = as<VectorXd>(tetav);
            MatrixXd kv = k_v[(*v-1)];
            norm_hv = (sqrt(tetv.transpose()* kv*tetv));

            Norm_h[(*v-1)] = norm_hv;
            fitted = fitted + fitv(_,(*v-1));
            penalite = penalite+(mu*norm_hv);
          }
        }
        Norm_H = Norm_h;
        NumericVector yfit; yfit = Y-fitted;
        SCR = pow(as<VectorXd>(yfit).norm(),2);
        crit = (SCR/sqrt(n))+penalite;
        rDiffCrit = (criti-crit)/criti;
        if(verbose){
          Rcout<<"Active support is "<<support<<"\n";
          Rcout<<"rDiffCrit=crit_i-crit_(i-1)/crit_i is "<<rDiffCrit<<"\n";
          Rcout<<"nrmRel=(norm(teta_i-teta_(i-1))/norm(teta_(i-1)))^2 is "<<nrmRel<<"\n";
        }
        // 0 is false and 1 is true : 0&1 = 0 ,...
        int brk; brk = ((nrmRel<eps)&(rDiffCrit<eps));
        // 0&&1 = false
        if (brk&&(i>=1)) break;
      }//if(convergence)
    } while (i<=maxIter);
    if (i>maxIter) {
      convergence = false;
    }
    colnames(fitv) = namesGrp;
    rownames(teta) = namesGrp;
    Norm_H.names() = namesGrp;
    StringVector Nsupp;
    if(support.size()>0) Nsupp = namesGrp[(support-1)];
    L = List::create(Named("intercept",c),Named("teta",teta)
                       ,Named("fit.v",fitv),Named("fitted",fitted)
                       ,Named("Norm.H",Norm_H)
                       ,Named("supp",support),Named("Nsupp",Nsupp)
                       ,Named("SCR",SCR),Named("crit",crit)
                       ,Named("MaxIter",i),Named("convergence",convergence)
                       , Named("RelDiffCrit",rDiffCrit),Named("RelDiffPar",nrmRel));
  }//End if(active.size()!=0)
  return L;
}//End RKHSgrplasso
