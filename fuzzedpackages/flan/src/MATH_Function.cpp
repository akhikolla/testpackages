
#include "MATH_Function.h"


using namespace Rcpp;


double MATH_Integration::computeFunction(double x, void* par){
// std::cout<<"Set the integrand"<<std::endl;
double res=0.;
double n0=0,n1=0,d0=0,d1=0;
if(mName.compare("CLONE_P0_WD") == 0){
  MATH_Params* params = (MATH_Params*) par;
  // std::cout<<"Set value of params"<<std::endl;
  double rho = params->rho;
  // std::cout<<"rho ="<<rho<<std::endl;
  double delta = params->delta;
  // std::cout<<"delta ="<<delta<<std::endl;
  // double k = params->k;
  res=(1-x)*pow(x,rho-1)/(1-delta*x);
}
else if(mName.compare("CLONE_P0_WD_WPEF") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  double zeta = params->zeta;

  n0=delta*(1-x)-(1-zeta)*(delta-x*(1-delta));
  d0=1-delta-delta*x-(1-zeta)*(1-delta)*(1-x);

  res=n0/d0*pow(x,rho-1);
}
else if(mName.compare("CLONE_PK_WD") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  double k = params->k;

  res=pow(x,rho)*pow(1-x,k-1)/pow(1-x*delta,k+1);
}
else if(mName.compare("CLONE_PK_WD_WPEF") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  double zeta = params->zeta;
  double k = params->k;

  n0=delta*(1-x)-(1-zeta)*(delta-x*(1-delta));
  n1=-zeta*(delta-x*(1-delta));
  d0=1-delta-delta*x-(1-zeta)*(1-delta)*(1-x);
  d1=-zeta*(1-delta)*(1-x);

  res=(n1*d0-n0*d1)/pow(d0,2.)*pow(x,rho-1);
  if(k > 1) res*=pow(-d1/d0,k-1);
}
else if(mName.compare("CLONE_dP0_dr_WD") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  // double k = params->k;

  res=(1-x)*pow(x,rho-1)/(1-delta*x)*log(x);

}
else if(mName.compare("CLONE_dP0_dr_WD_WPEF") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  double zeta = params->zeta;

  n0=delta*(1-x)-(1-zeta)*(delta-x*(1-delta));
  d0=1-delta-delta*x-(1-zeta)*(1-delta)*(1-x);

  res=n0/d0*pow(x,rho-1)*log(x);
}
else if(mName.compare("CLONE_dPK_dr_WD") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  double k = params->k;

  res=pow(x,rho)*pow(1-x,k-1)/pow(1-x*delta,k+1)*log(x);
}
else if(mName.compare("CLONE_dPK_dr_WD_WPEF") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  double zeta = params->zeta;
  double k = params->k;

  n0=delta*(1-x)-(1-zeta)*(delta-x*(1-delta));
  n1=-zeta*(delta-x*(1-delta));
  d0=1-delta-delta*x-(1-zeta)*(1-delta)*(1-x);
  d1=-zeta*(1-delta)*(1-x);

  res=(n1*d0-n0*d1)/pow(d0,2.)*pow(x,rho-1)*log(x);
  if(k > 1) res*=pow(-d1/d0,k-1);
}
else if (mName.compare("CLONE_PGF") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  // double k = params->k;

  res=pow(x,rho)/(1+x*delta);
}
else if (mName.compare("CLONE_dPGF_dr") == 0){
  MATH_Params* params = (MATH_Params*) par;
  double rho = params->rho;
  double delta = params->delta;
  // double k = params->k;

  res=pow(x,rho)/(1+x*delta)*log(x);
}
// else if(mName.compare("TEST") == 0){
//   MATH_Params* params = (MATH_Params*) par;
//   double rho = params->rho;
//   // res=pow(x,rho)*log(x);
//   res=pow(x,rho)/pow(1-x*delta,2)*log(x);
//
// }
// std::cout<<"Set the parameters"<<std::endl;
  return(res);
}

double gslClassWrapper(double x, void * pp) {
    MATH_Params *p = (MATH_Params *)pp;
    return p->ptMInt->computeFunction(x,p);
}


double MATH_Integration::computeIntegral(double a, double b) {

//   List res;
//   if(k > 0) res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["subdivisions"] = mSubd,_["rho"] =rho,_["delta"] =delta,_["k"]=k);
//   else res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["subdivisions"] = mSubd,_["rho"] =rho,_["delta"] =delta);
//   double integ=as<double>(res["value"]);
// //   double integ=res["value"];
//   return  integ;
mParams->ptMInt = this;

gsl_function integrand;

integrand.params=mParams;
integrand.function = &gslClassWrapper;

// std::cout<<"Reserve workspace"<<std::endl;
gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
double result, error;
// std::cout<<"Compute integral"<<std::endl;
gsl_integration_qags(&integrand, a, b, 0, mReltol, mSubd,w, &result, &error);
// std::cout<<"Free workspace"<<std::endl;
gsl_integration_workspace_free(w);

return result;

}



void MATH_Polynom::reduce(double eps){
      int i=0,dmax=mCoef.size()-1;
      std::vector<double>::iterator it1=mCoef.begin();
      std::vector<double>::iterator itmax;
      for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it,i++) {
	if(*it <= eps) *it=0;
	if(*it > 0) {
	  itmax=it;
	  dmax=i ;
	}

      }

      mCoef=std::vector<double>(it1,itmax+1);
      mDeg=dmax;
    }

// #define REAL(z,i) ((z)[2*(i)])
// #define IMAG(z,i) ((z)[2*(i)+1])

void MATH_Polynom::square_fft(){

  int n=mCoef.size();

  //RcppGSL version
  // gsl_fft_complex_wavetable * wavetable;
  // gsl_fft_complex_workspace * workspace;
  // // std::cout<<"Initial size of mCoef : n="<<n<<std::endl;
  //
  // std::vector<double> temp=mCoef;
  // n*=2;
  // n--;
  //   // std::cout<<"New value of n ="<<n<<std::endl;
  // temp.resize(n);
  // double tp[2*n];
  // gsl_complex_packed_array cx_coef=tp;
  // for(int i=0;i<n;i++){
  //   REAL(cx_coef,i)=temp[i];
  //   IMAG(cx_coef,i)=0.;
  // }
  //
  // wavetable = gsl_fft_complex_wavetable_alloc (n);
  // workspace = gsl_fft_complex_workspace_alloc (n);
  //
  // // FFT
  // std::cout<<"FFT"<<std::endl;
  // gsl_fft_complex_forward (cx_coef, 1, n,wavetable, workspace);
  // // temp=cx_coef;
  // double r,im;
  // // std::cout<<"Carré des coeff après fft"<<std::endl;
  //
  // for(int i=0; i<n;i++){
  //   r=REAL(cx_coef,i);
  //   im=IMAG(cx_coef,i);
  //   REAL(cx_coef,i) = r*r-im*im;
  //   IMAG(cx_coef,i) = 2*r*im;
  // }
  // // IFFT
  // std::cout<<"FFT Inverse"<<std::endl;
  // gsl_fft_complex_inverse (cx_coef, 1, n,wavetable, workspace);
  //
  // mCoef.resize(n);
  // mDeg=n-1;
  //
  // for(int i=0; i<n;i++){
  //   mCoef[i]=REAL(cx_coef,i);
  // }
  // gsl_fft_complex_wavetable_free (wavetable);
  // gsl_fft_complex_workspace_free (workspace);

  //RcppArmadillo version
  std::vector<double> temp=mCoef;
  n*=2;
  n--;
  // std::cout<<"New value of n ="<<n<<std::endl;
  temp.resize(n);

  double tp=log2((double)(n));


  int N=floor(tp);

  // std::cout<<"Initial value of N="<<N<<std::endl;


  if(N != tp) {
    N++;
    // std::cout<<"New value of N ="<<N<<std::endl;
    N=pow(2.,N);
    // std::cout<<"Final value of N ="<<N<<std::endl;
    temp.resize(N);
  }
  // std::cout<<"Final value of N ="<<N<<std::endl;

//   ComplexVector CFFT=(*mFFT)(temp,_["inverse"]=false);
  arma::vec tempp(temp);
  arma::cx_vec CFFT=arma::fft(tempp);

  arma::cx_double c(0,0);
  arma::cx_vec::iterator it;
  // std::cout<<"Carré des coeff après fft"<<std::endl;
  for(it=CFFT.begin() ; it!=CFFT.end(); ++it) {
    // std::cout<<*it<<std::endl;
    c=*it;
    // std::cout<<"devient"<<std::endl;
    *it=c*c;
    // std::cout<<*it<<std::endl;
  }

  CFFT=arma::ifft(CFFT);
//   CFFT = (*mFFT)(CFFT,_["inverse"]=true);

  mCoef.resize(n);
  mDeg=n-1;

  it=CFFT.begin();
  // std::cout<<"Coeff après ifft"<<std::endl;
  for(std::vector<double>::iterator itC=mCoef.begin() ; itC!=mCoef.end() ; ++itC , ++it){
//     *itC = (*it).r;
    *itC = std::real(*it);
    // std::cout<<*itC<<std::endl;
    // (*itC)/=N;
  }
}
