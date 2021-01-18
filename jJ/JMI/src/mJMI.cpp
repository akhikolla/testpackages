#include <vector>
#include <map>
#include <algorithm>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;
arma::mat tiedrank(arma::mat x){
  int nrow=x.n_rows,ncol=x.n_cols;
  arma::mat x_sorted=sort(x);
  arma::mat t_rank;t_rank.zeros(nrow,ncol);
  arma::mat rank;
  rank.zeros(x.n_rows,x.n_cols);
  for(int i=0;i<ncol;i++){

    map<double, double> val2t_rank;
    double flag=x_sorted(0,i);
    int ppos=0;
    for(int j=1;j<nrow;j++){
      if((abs(x_sorted(j,i)-flag)>1e-12)){
        for(int k=ppos;k<j;k++){
          rank(k,i)=(ppos+j-1.0)/2;
          val2t_rank.insert(make_pair(x_sorted(k,i),(ppos+j-1.0)/2));
        }
        flag=x_sorted(j,i);
        ppos=j;
      }
    }
    int j = x.n_rows;
    for(int k=ppos;k<j;k++){
      rank(k,i)=(ppos+j-1.0)/2;
      val2t_rank.insert(make_pair(x_sorted(k,i),(ppos+j-1.0)/2));
    }

    for(int j=0;j<nrow;j++){
      t_rank(j,i)=val2t_rank[x(j,i)];
    }
  }

  return t_rank;
}



//[[Rcpp::export]]
List mJMICpp(const arma::mat& x,const arma::mat& y, const int BN){
  int n = x.n_rows, p=x.n_cols;
  int q = y.n_cols;


  arma::umat BS(n,BN);
  for (int ip=0;ip<BN;ip++){
    vec ss=randu<vec>(n);
    uvec ss1=sort_index(ss);
    BS.col(ip)=ss1;    
  }

  arma::mat nx = tiedrank(x)/(n+1);
  arma::mat ny = tiedrank(y)/(n+1);
  arma::cube kernelx2;
  kernelx2.zeros(n,n,p);
  for(int ip=0;ip<p;ip++){
    arma::mat kx=repmat(nx.col(ip),1,n);
    kernelx2.slice(ip)=abs(kx-trans(kx));
  }
  arma::cube kernelx21=kernelx2;
  arma::cube kernelx22=pow(kernelx2,2.0);
  arma::cube kernelx23=pow(kernelx2,3.0);
  arma::cube kernelx24=pow(kernelx2,4.0);

  arma::cube kernely2;
  kernely2.zeros(n,n,q);
  for(int iq=0;iq<q;iq++){
    arma::mat ky=repmat(ny.col(iq),1,n);
    kernely2.slice(iq)=abs(ky-trans(ky));
  }
  arma::cube kernely21=kernely2;
  arma::cube kernely22=pow(kernely2,2.0);
  arma::cube kernely23=pow(kernely2,3.0);
  arma::cube kernely24=pow(kernely2,4.0);

  double sigma=sqrt(n/(12*(n+1.0)));
  double h=sigma/(pow(n,1/(p+q+3.0)));
  static const double H[50] = { 0.02, 0.08, 0.18, 0.32, 0.50000000000000011,
                                  0.72, 0.98000000000000009, 1.28, 1.6199999999999999, 2.0000000000000004,
                                  2.42, 2.88, 3.3800000000000003, 3.9200000000000004, 4.5, 5.12,
                                  5.7800000000000011, 6.4799999999999995, 7.22, 8.0000000000000018,
                                  8.8199999999999985, 9.68, 10.58, 11.52, 12.5, 13.520000000000001,
                                  14.580000000000002, 15.680000000000001, 16.82, 18.0, 19.220000000000002,
                                  20.48, 21.78, 23.120000000000005, 24.499999999999996, 25.919999999999998,
                                  27.38, 28.88, 30.42, 32.000000000000007, 33.62, 35.279999999999994, 36.98,
                                  38.72, 40.5, 42.32, 44.18, 46.08, 48.019999999999996, 50.0 };
  int nH=50;
  arma::colvec mi(50);
  arma::colvec Fi(50);
  arma::mat BMI;
  BMI.zeros(BN,50);
  for(int iw=0;iw<nH;iw++){
  double bw=H[iw];
  double hw=bw*h;
  arma::cube kernelx = 1/(1.0+4*kernelx21/hw+6*kernelx22/(pow(hw,2))+4*kernelx23/(pow(hw,3))+kernelx24/(pow(hw,4)));
  arma::mat kx;
  kx.ones(n,n);
  for (int ip=0; ip<p;ip++){
    kx=kernelx.slice(ip)%kx;
  }
  arma::mat kx1=kx;
  for (int i=0;i<n;i++){
    kx(i,i)=0;
  }
  arma::colvec fx=mean(kx,1);
  arma::colvec fx1=mean(kx1,1);

  arma::cube kernely = 1/(1.0+4.0*kernely21/hw+6.0*kernely22/(pow(hw,2.0))+4.0*kernely23/(pow(hw,3.0))+kernely24/(pow(hw,4.0)));
  arma::mat ky;
  ky.ones(n,n);
  for (int iq=0; iq<q;iq++){
    ky=kernely.slice(iq)%ky;
  }
  arma::mat ky1=ky;
  for (int i=0;i<n;i++){
    ky(i,i)=0;
  }
  arma::colvec fy=mean(ky,1);
  arma::colvec fy1=mean(ky1,1);

  arma::colvec fxfy=fx%fy+1e-100;
  arma::colvec fxfy1=fx1%fy1+1e-100;
  arma::colvec A=log(mean(kx%ky,1)/fxfy);
  arma::colvec A1=log(mean(kx1%ky1,1)/fxfy1);
  mi[iw]=mean(A);
  Fi[iw]=mean(A1);

  for (int j=0;j<BN;j++){
  arma::mat kyj = ky(BS.col(j),BS.col(j));
  arma::colvec fyj=mean(kyj,1);

  arma::colvec fxfyj=fx%fyj+1e-100;
  arma::colvec Aj=log(mean(kx%kyj,1)/fxfyj);
    BMI(j,iw) = mean(Aj);
  }
  }

  colvec B = max(BMI,1);
  int locate=index_max(mi);
  double c=0.5*Fi(locate);
  double jmila=mi[locate];
  mi=0.5*jmila+c;
  int counter=0;double pvalue=0;
  for (int sss=0;sss<BN;sss++){
    if (B[sss]>jmila){
      counter++;
    }
  }
  pvalue=(counter+0.0)/BN;
  double jmi=max(mi);
  //return jmif,palue
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("mi")=jmi,
                                  Rcpp::Named("pvalue")=pvalue);
  return result;
}
