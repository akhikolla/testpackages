#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <math.h>
#include <gsl/gsl_integration.h>
// #include <gsl/gsl_fft_complex.h>


using namespace Rcpp;


// forward definition of Parameters struct
struct MATH_Params;

// Class definition of Integrator
class MATH_Integration {

private:

  double mReltol;

  int mSubd;

  // Function* mIntegrate=NULL;

  // Function* mIntegrand=NULL;
  MATH_Params* mParams=NULL;
  std::string mName;

  // List mFcts;


protected:

public:

    MATH_Integration() {
      // if(!mParams) delete mParams;
      // mParams=NULL;
    };
    // MATH_Integration(List fns, double reltol, int subd){
    MATH_Integration(double reltol, int subd){
      mReltol=reltol;
      mSubd=subd;
      //
      // if(!mParams) delete mParams;
      // mParams=NULL;
    };


    ~MATH_Integration(){
      // if(!mParams) delete mParams;
      // if(!mIntegrate) delete mIntegrate;
      // if(!mIntegrand) delete mIntegrand;
    };

    double computeFunction(double x, void* par);
    void setFunctionName(std::string name) {
      mName=name;
    };
    void setFunction(std::string name,MATH_Params* par){
      // double rho,double delta, double k){
      mName=name;
      mParams=par;
    };

  /*
   * Integrals functions
   */

  double computeIntegral(double a, double b);

};

// Parameters struct
struct MATH_Params {
  double rho;
  double delta;
  double zeta;
  double k;
  MATH_Integration* ptMInt;

};



// double gslClassWrapper(double x, void * pp) ;


class MATH_Polynom {

private:

  void init_all(int deg){
    setDegree(deg);
//     init_fft();
  };

//   void init_fft(){
//     if(!mFFT) delete mFFT;
//     mFFT= new Function("fft");
//   };

protected:
   std::vector<double> mCoef;
   int mDeg;

//    Function* mFFT=NULL;


public:

    MATH_Polynom() {
      init_all(1);
    };

    MATH_Polynom(std::vector<double> C){
      setCoef(C);
//       init_fft();
    }

    ~MATH_Polynom(){
//       if(!mFFT) delete mFFT;
    };


    void setDegree(int deg){
      mCoef.resize(deg+1);
      mCoef[deg]=1;
      mDeg=deg;
    };

    void setCoef(std::vector<double> C){
      mCoef=C;
      mDeg=C.size()-1;
    };


    std::vector<double> getCoef(){
      return mCoef;
    }

    int getDegree(){
      return mDeg;
    };

    double& operator[](int i) {
      return mCoef[i];
    };

    const double& operator[](int i) const {
      return mCoef[i];
    };

    MATH_Polynom& operator*=(double f) {
      int i=0;
      for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it,i++) {
	(*it)*=f;
      }
      return *this;
    };

    MATH_Polynom& operator+=(double f) {
      mCoef[0]+=f;
      return *this;
    };

    void reduce(double eps);

    void square_fft();


};



#endif
