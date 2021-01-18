#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace Rcpp;
NumericVector k0linearT(NumericVector u, double v) {
  int l; l = u.length();
  NumericVector k(l);
  NumericVector k0(l);
  k = 1+v*u;
  k0 = k - 4*(1+v/2)*(1+u/2)/5;
  return(k0);
}//End k0linearT
NumericVector k0quadT(NumericVector u, double v) {
  int l; l = u.length();
  Rcpp::NumericVector k(l);
  Rcpp::NumericVector k0(l);
  Rcpp::NumericVector uv(l);
  uv = 1+(v*u);
  k = pow(uv,2);
  Rcpp::NumericVector uu(l);
  double vv;
  uu = pow(u,2)/3;
  vv = v*v/3;
  double a;
  Rcpp::NumericVector b(l);
  Rcpp::NumericVector c(l);
  a=1+v+vv;
  b=1+u+uu;
  Rcpp::NumericVector ab(l);
  ab=a*b;
  c=3*ab/5;
  k0 = k - c;
  return(k0);
}//End k0quadT
NumericVector zmT(NumericVector x){
  int l; l = x.length();
  Rcpp::NumericVector z(l);
  z = 2-(1+x)*exp(-2*x)-(2-x)*exp(-2*(1-x));
  return (z);
}
double zmdT(double x){
  double zz;
  zz = 2-(1+x)*exp(-2*x)-(2-x)*exp(-2*(1-x));
  return (zz);
}
NumericVector k0maternT(NumericVector u, double v) {
  int l; l = u.length();
  Rcpp::NumericVector k(l);
  Rcpp::NumericVector k0(l);
  k = (1+2*abs(u-v))*exp(-2*abs(u-v));
  const double z2 = 0.8383385;
  k0 = k - zmT(u)*zmdT(v)/z2;
  return(k0);
}//End k0maternT
NumericVector k0brownianT(NumericVector u, double v) {
  int l; l = u.length();
  NumericVector k(l);
  NumericVector k0(l);
  double vb;
  k = 1+pmin(u,rep(v,u.length()));
  vb=1+v*(1-v/2);
  k0 = k-3*(1+u*(1-u/2))*vb/4;
  return(k0);
}//End k0brownianT
NumericVector zgT(NumericVector x){
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
double zgdT(double x){
  double p1d;
  double p2d;
  double zgd;
  p1d = R::pnorm((1-x)*1.414214,0,1,TRUE,FALSE);
  p2d =  R::pnorm(-x*1.414214,0,1,TRUE,FALSE);
  zgd = 1.772454*(p1d-p2d);
  return (zgd);
}
NumericVector k0gaussianT(NumericVector u, double v) {
  int n;
  n = u.size();
  Rcpp::NumericVector k(n);
  Rcpp::NumericVector k0(n);
  Rcpp::NumericVector p(n);
  p = pow((u-v),2);
  k = exp(-p);
  const double z2g = 0.8615277;
  k0 = k - zgT(u)*zgdT(v)/z2g;
  return(k0);
}//End k0gaussianT
NumericVector k0T(NumericVector u, double v, String kernel) {
  int l; l = u.length();// I add length
  Rcpp::NumericVector k(l);
  if (kernel=="matern") k = k0maternT(u,v);
  if (kernel=="brownian") k = k0brownianT(u,v);
  if (kernel=="gaussian") k = k0gaussianT(u,v);
  if (kernel=="linear") k = k0linearT(u,v);
  if (kernel=="quad") k = k0quadT(u,v);
  return(k);
}//End k0T
SEXP KvTest(NumericMatrix X,NumericMatrix XT,String kernel,int Dmax){
  List index_K_T(Dmax);
  int d = X.cols();
  int n = X.rows();
  int nT = XT.rows();
  unsigned int nnT=n*nT;
  int lz=0;
  for(int suz=1; suz<Dmax+1 ; suz++){
    lz+=Rf_choose(d,suz);
  }
  List matZ(lz);

  for (int hh=0; hh<d; hh++) {
    int hk=0;
    NumericVector k(nnT);
    for (int jj=0; jj<n; jj++) {
      NumericVector v0k(nT);
      v0k = k0T(XT(_,hh),X(jj,hh),kernel);
      std::copy(v0k.begin(), v0k.end() ,k.begin()+(nT*hk));
      hk+=1;
    }
    k.attr("dim") = Dimension(nT, n);
    NumericMatrix mk; mk = as<NumericMatrix>(k);
    matZ[hh]=mk;
  }
  NumericVector index0(d);
  index0=seq(1,d);
  index0.attr("dim")=Dimension(d,1);
  index_K_T[0] = index0;
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
        NumericVector vj;
        NumericVector vl;
        vj = matZ[j-1];
        vl = matZ[l-1];
        NumericVector vjvl(nnT);
        vjvl = vj*vl;
        index2_v1[r]=j;
        index2_v2[r]=l;
        r+=1;
        vjvl.attr("dim") = Dimension(nT, n);
        NumericMatrix mv; mv = as<NumericMatrix>(vjvl);
        matZ[d+h]=mv;
        h+=1;
      }
    }
    NumericMatrix index2(nrowindex2,2);
    index2 = cbind(index2_v1,index2_v2);
    index_K_T[i] = index2;
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
              NumericVector vjj;
              NumericVector vhh;
              int indexjj = jj+clk-c-1;
              vjj = matZ[indexjj];
              vhh = matZ[hh];
              NumericVector vjjvhh(nnT);
              vjjvhh = vjj*vhh;
              vjjvhh.attr("dim") = Dimension(nT, n);
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
    }
  }
  return matZ;
}//End KvTest
// [[Rcpp::export]]
NumericMatrix PredErr(NumericMatrix X, NumericMatrix XT, NumericVector YT,NumericVector mu,
                              NumericVector gamma, List res, String kernel, int Dmax){
  List reskv; reskv = KvTest(X,XT,kernel,Dmax);
  int l; l = res.size();
  int nT; nT = XT.rows();
  int n; n = X.rows();
  int lm; lm = mu.size();
  int lg; lg = gamma.size();
  NumericVector verr(lm*lg);
    for(int ll=0; ll<l;ll++){
      List pen;
      pen = res[ll];

      List metmod; metmod = pen[2];
      double intercept; intercept = metmod["intercept"];
      NumericVector fitted(nT); fitted = rep(intercept,nT);
      NumericVector supp; supp = metmod["supp"];
      int ls; ls = supp.size();
      if(ls>0){
        MatrixXd teta; teta = metmod["teta"];
        for(NumericVector::iterator v = supp.begin(); v != supp.end(); ++v){
          VectorXd tetav(n); tetav = teta.row(*v-1);
          MatrixXd kv(nT,n); kv = reskv[*v-1];
          NumericVector kvtetav(nT); kvtetav = kv*tetav;
          fitted = fitted + kvtetav;
        }
      }//if(ls>0)
      verr(ll) = sum(pow((YT-fitted),2))/nT;
    }//for(int ll=0; ll<l;ll++)
  verr.attr("dim") = Dimension(lg,lm);
    NumericMatrix ErrPred; ErrPred = as<NumericMatrix>(verr);
    StringVector a;
    for(int i=1;i<lm+1;i++){
      String b("mu = ");
      a.push_back(b+=mu[i-1]);
    }
    StringVector aa;
    for(int j=1;j<lg+1;j++){
      String bb("gamma = ");
      aa.push_back(bb+=gamma[j-1]);
    }
    colnames(ErrPred)=a;
    rownames(ErrPred)=aa;
  return ErrPred;
}//End PredErr
