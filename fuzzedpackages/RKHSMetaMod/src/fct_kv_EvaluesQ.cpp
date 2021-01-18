#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
#include "fct_kv_EvaluesQ.h"
NumericVector k0linear(NumericVector u, double v) {
  int l; l = u.length();
  NumericVector k(l);
  NumericVector k0(l);
  k = 1+v*u;
  k0 = k - 4*(1+v/2)*(1+u/2)/5;
  return(k0);
}//End k0linear
NumericVector k0quad(NumericVector u, double v) {
  int l; l = u.length();
  NumericVector k(l);
  NumericVector k0(l);
  NumericVector uv(l);
  uv = 1+(v*u);
  k = pow(uv,2);
  NumericVector uu(l);
  double vv;
  uu = pow(u,2)/3;
  vv = v*v/3;
  double a;
  NumericVector b(l);
  NumericVector c(l);
  a=1+v+vv;
  b=1+u+uu;
  NumericVector ab(l);
  ab=a*b;
  c=3*ab/5;
  k0 = k - c;
  return(k0);
}//End k0quad
NumericVector int_1v(NumericVector x){
  int l; l = x.length();
  NumericVector z(l);
  z = 2-(1+x)*exp(-2*x)-(2-x)*exp(-2*(1-x));
  return (z);
}
double int_1d(double x){
  double zz;
  zz = 2-(1+x)*exp(-2*x)-(2-x)*exp(-2*(1-x));
  return (zz);
}
NumericVector k0matern(NumericVector u, double v) {
  int l; l = u.length();
  NumericVector k(l);
  NumericVector k0(l);
  k = (1+2*abs(u-v))*exp(-2*abs(u-v));
  const double int_2 = 0.8383382;
  k0 = k - int_1v(u)*int_1d(v)/int_2;
  return(k0);
}//End k0matern
NumericVector k0brownian(NumericVector u, double v) {
  int l; l = u.length();
  NumericVector k(l);
  NumericVector k0(l);
  double vb;
  k = 1+pmin(u,rep(v,u.length()));
  vb=1+v*(1-v/2);
  k0 = k-3*(1+u*(1-u/2))*vb/4;
  return(k0);
}//End k0brownian
NumericVector int_1gv(NumericVector x){
  int n;
  n = x.size();
  NumericVector p1v(n);
  NumericVector p2v(n);
  NumericVector zgv(n);
  p1v = pnorm((1-x)*1.414214);
  p2v =  pnorm(-x*1.414214);
  zgv = 1.772454*(p1v-p2v);
  return (zgv);
}
double int_1gd(double x){
  double p1d;
  double p2d;
  double zgd;
  p1d = R::pnorm((1-x)*1.414214,0,1,TRUE,FALSE);
  p2d = R::pnorm(-x*1.414214,0,1,TRUE,FALSE);
  zgd = 1.772454*(p1d-p2d);
  return (zgd);
}
NumericVector k0gaussian(NumericVector u, double v) {
  int n;
  n = u.size();
  NumericVector k(n);
  NumericVector k0(n);
  NumericVector p(n);
  p = pow((u-v),2);
  k = exp(-p);
  const double int_2g = 0.8615277;
  k0 = k - int_1gv(u)*int_1gd(v)/int_2g;
  return(k0);
}//End k0gaussian
NumericVector k0(NumericVector u, double v, String kernel) {
  int l; l = u.length();// I add length
  Rcpp::NumericVector k(l);
  if (kernel=="matern") k = k0matern(u,v);
  if (kernel=="brownian") k = k0brownian(u,v);
  if (kernel=="gaussian") k = k0gaussian(u,v);
  if (kernel=="linear") k = k0linear(u,v);
  if (kernel=="quad") k = k0quad(u,v);
  return(k);
}//End k0
StringVector concatenate(StringVector a){
  StringVector c;
  std::ostringstream x;
  for (int i = 0; i < a.size(); i++)
  x << a[i];
  c.push_back(x.str());
  return c;
}//End concatenate
StringVector namesGrp(int d, int Dmax, List index){
  StringVector a;
  for(int i=1;i<d+1;i++){
    String b("v");
    a.push_back(b+=i);
  }
  for(int ij=0;ij<d;ij++){
    a[ij]+=".";
  }
  if (Dmax >1) {
    for (int D=2;D<Dmax+1;D++) {
      NumericVector setD;
      setD = index[D-1];
      NumericVector dimsetD(2);
      dimsetD(0)=Rf_choose(d,D) ; dimsetD(1)= D;
      setD.attr("dim")= dimsetD;
      NumericMatrix m;
      m = as<NumericMatrix>(setD);

      for (int ic=0;ic<Rf_choose(d,D);ic++) {
        String b("v");
        IntegerVector vv;
        vv = m(ic,_);
        CharacterVector y = as<CharacterVector>(vv);
        int lvv; lvv=vv.length();
        for(int ivv=0;ivv<lvv;ivv++){
          y[ivv]+=".";
        }
        StringVector bbnew(1); bbnew=concatenate(y);
        String bnew; bnew=b+bbnew[0];
        a.push_back(bnew);
      }
    }
  }
  return(a);
}//End namesGrp
// [[Rcpp::export]]
SEXP calc_Kv(NumericMatrix X, String kernel, int Dmax, bool correction=true, 
             bool verbose=true, double tol=1e-8){
  List L;
  List index_K_T(Dmax);
  StringVector namGrp;
  int d; d = X.cols();
  int n; n = X.rows();
  unsigned int nn=n*n;
  int lz=0;
  for(int suz=1; suz<Dmax+1 ; suz++){
    lz+=Rf_choose(d,suz);
  }
  List matZ(lz);
  for (int hh=0; hh<d; hh++) {
    int hk=0;
    NumericVector k(nn);
    for (int jj=0; jj<n; jj++) {
      NumericVector v0k(n);
      v0k = k0(X(_,hh),X(jj,hh),kernel);
      std::copy(v0k.begin(), v0k.end() ,k.begin()+(n*hk));
      hk+=1;
    }
    k.attr("dim") = Dimension(n, n);
    NumericMatrix mk; mk = as<NumericMatrix>(k);
    matZ[hh]=mk;
  }
  NumericVector index0(d);
  index0=seq(1,d);
  index0.attr("dim")=Dimension(d,1);
  index_K_T[0] = index0;

  /*case Dmax=1*/
  if(Dmax==1){
    namGrp=namesGrp(d,Dmax,index_K_T);
    matZ.attr("names")=namGrp;
  }

  /*case Dmax!=1*/
  if(Dmax!=1){
    int i = 1;
    int nrowindex2 = Rf_choose(d,2);
    NumericVector index2_v1(nrowindex2);
    NumericVector index2_v2(nrowindex2);
    int h=0;
    int r=0;
    for(int j=1;j<d;++j){
      for(int l=j+1;l<d+1;++l){
        NumericVector vj(n);
        NumericVector vl(n);
        vj = matZ[j-1];
        vl = matZ[l-1];
        NumericVector vjvl(nn);
        vjvl = vj*vl;
        index2_v1[r]=j;
        index2_v2[r]=l;
        r+=1;
        vjvl.attr("dim") = Dimension(n, n);
        NumericMatrix mv; mv = as<NumericMatrix>(vjvl);
        matZ[d+h]=mv;
        h+=1;
      }
    }
    NumericMatrix index2(nrowindex2,2);
    index2 = cbind(index2_v1,index2_v2);
    index_K_T[i] = index2;

    /* case Dmax=2 */
    if(Dmax==2){
      namGrp=namesGrp(d,Dmax,index_K_T);
      matZ.attr("names")=namGrp;
    }

    /* case Dmax>2 */
    if(Dmax>2){
      for(int ii=2;ii<Dmax;++ii){
        NumericMatrix mm_index(Rf_choose(d,(ii+1)),(ii+1));
        int st=0;
        NumericMatrix m_dd=index_K_T[ii-1];
        int dim=ii;
        int c=Rf_choose(d,ii);
        int clk=0;
        for(int csu=1; csu<ii+1 ; csu++){
          clk+=Rf_choose(d,csu);
        }
        int s=0;
        for(int jj=1;jj<(c+1);++jj){
          int dd;
          dd=m_dd(jj-1,dim-1);
          if(dd<d){
            for(int hh=dd;hh<d;++hh){
              NumericVector vjj(n);
              NumericVector vhh(n);
              int indexjj = jj+clk-c-1;
              vjj = matZ[indexjj];
              vhh = matZ[hh];
              NumericVector vjjvhh(nn);
              vjjvhh = vjj*vhh;
              vjjvhh.attr("dim") = Dimension(n, n);
              NumericMatrix mvi; mvi = as<NumericMatrix>(vjjvhh);
              matZ[clk+s]=mvi;
              s+=1;
            }
          }
          int ss=0;
          int index=d-dd;
          if(index!=0){
            NumericVector index1=m_dd((jj-1),_);
            NumericVector index22(d-dd);
            index22=seq((dd+1),d);
            int leng=index22.length();
            NumericVector index11=rep_each(index1,leng);
            for(int l=0;l<(ii+1);++l){
              if(l!=ii){
                for(int ll=0;ll<index22.length();++ll){
                  mm_index(st+ll,l)=index11[ss];
                  ss+=1;
                }
              }
              if(l==ii){
                ss=0;
                for(int ll=0;ll<index22.length();++ll){
                  mm_index(st+ll,l)=index22[ss];
                  ss+=1;
                }
              }
            }
            st+=leng;
          }
        }
        index_K_T[ii]=mm_index;
      }
      namGrp=namesGrp(d,Dmax,index_K_T);
      matZ.attr("names")=namGrp;
    }
  }
  // Second step : calculate the eigenvalues and eigenvectors of kvs
    for(int ij=0;ij<lz;ij++){
      List r(3);
      MatrixXd Mkv; Mkv=matZ[ij];
      SelfAdjointEigenSolver<MatrixXd> es(Mkv);
      VectorXd eV = es.eigenvalues();
      MatrixXd V = es.eigenvectors();
      if(correction){
        double dmax; dmax = eV.maxCoeff();
        double epsilonk; epsilonk = tol*dmax;
        double dmin; dmin = eV.minCoeff();
        if(dmin <= epsilonk){
        for(int i=0;i<eV.size();i++){
          eV[i]+=epsilonk;
        }
        StringVector nv; nv=namGrp;
        if(verbose){
          Rcout<<"The nearest positive definite matrix is calculated for group :"<<nv[ij]<<"\n";
        }
      }
    }//if(correction)
      r = List::create(Named("Evalues",eV),Named("Q",V));
      matZ[ij]=r;
  }//for(int ij=0;ij<lz;ij++)
  L = List::create(Named("kv",matZ),Named("names.Grp",namGrp));
  return L;
}//End calc_Kv
