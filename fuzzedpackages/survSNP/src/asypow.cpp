#include "integration_wrapper.h"
#include <Rcpp.h>

template<class K> 
class FunctionMaker {
public:
    static K* a;
    static double (K::*func)(double);
    static double funcForIntegration(double x) {
        return (a->*func)(x);
    }
};
template<class K> K* FunctionMaker<K>::a = NULL;
template<class K> double (K::*FunctionMaker<K>::func)(double) = NULL;

double mean(double * a, int n) {
    double sum=0;
    for(int i=0;i<n;i++) {
         sum+=a[i];
    }
    return sum/((double)n);
}

double pnormStd(double x) {
     return Rcpp::as<double>(Rcpp::wrap(  Rcpp::pnorm(Rcpp::NumericVector::create(x),0.0, 1.0)  ));
}

double qnormStd(double p) {
     return Rcpp::as<double>(Rcpp::wrap(  Rcpp::qnorm(Rcpp::NumericVector::create(p),0.0, 1.0)  ));
}

class Asypow {
public:
    double n, theta,a,b,lambda0,q,alpha;
    Rcpp::NumericVector p,z;
    bool exactvar;

    Asypow(double _n, double _theta, double _a, double _b, double _lambda0, double _q, 
           Rcpp::NumericVector _p, double _alpha, Rcpp::NumericVector _z, bool _exactvar)
    : n(_n), theta(_theta), a(_a), b(_b), lambda0(_lambda0), q(_q), alpha(_alpha), p(_p), z(_z), exactvar(_exactvar)
    {}

    void setValue(double _n, double _theta, double _a, double _b, double _lambda0, double _q,
                  Rcpp::NumericVector _p, double _alpha, Rcpp::NumericVector _z, bool _exactvar)
    {
        n=_n;
        theta=_theta;
        a=_a;
        b=_b;
        lambda0=_lambda0;
        q=_q;
        p=_p;
        alpha=_alpha;
        z=_z;
        exactvar=_exactvar;
    }
    ///////////////////////////////////////////////////////
    double s_c(double t){
        return (t<a)+(t>=a&&t<=b)*(1-(t-a)/(b-a));
    }
    double s0(double t){
        return ((1-q)*(1-q)*exp(-lambda0*t)
               +2*q*(1-q)*exp(theta)*exp(-exp(theta)*lambda0*t)
               +q*q*exp(2*theta)*exp(-exp(2*theta)*lambda0*t))*s_c(t);
    }
    double s1(double t){
        return (2*q*(1-q)*exp(theta)*exp(-exp(theta)*lambda0*t)
                +2*q*q*exp(2*theta)*exp(-exp(2*theta)*lambda0*t))*s_c(t);
    }
    double s2(double t){
        return (2*q*(1-q)*exp(theta)*exp(-exp(theta)*lambda0*t)
               +4*q*q*exp(2*theta)*exp(-exp(2*theta)*lambda0*t))*s_c(t);
    }
    double r0(double t){
        return ((1-q)*(1-q)*exp(-lambda0*t)
               +2*q*(1-q)*exp(-exp(theta)*lambda0*t)
               +q*q*exp(-exp(2*theta)*lambda0*t))*s_c(t);
    }
    double r1(double t){
        return (2*q*(1-q)*exp(-exp(theta)*lambda0*t)
               +2*q*q*exp(-exp(2*theta)*lambda0*t))*s_c(t);
    }
    double r2(double t){
        return (2*q*(1-q)*exp(-exp(theta)*lambda0*t)
               +4*q*q*exp(-exp(2*theta)*lambda0*t))*s_c(t);
    }
    double e0(double t){
        return r1(t)/r0(t);
    }
    double e(double t){
        return s1(t)/s0(t);
    }
    double a1(double t){
        return -e0(t);
    }
    double a2(double t){
         return -s0(t)/r0(t);
    }
    double a3(double t){
         return s0(t)*r1(t)/(r0(t)*r0(t));
    }
    // asymptotic mean
    double f(double t){
        return e0(t)*s0(t);
    }
    //calculate I_0(\theta)
    double g(double s){
        return (r2(s)/r0(s)-(r1(s)/r0(s))*(r1(s)/r0(s)))*s0(s)*lambda0;
    }
    // asymptotic variance
    // variance of xi(theta)
    template<int N>
    double fN(double s){
        return (z[N]-e0(s))*(z[N]-e0(s))*exp(-exp(theta*z[N])*lambda0*s)*s_c(s);
    }
    template<int N>
    double vN(){
        return p[N]*exp(theta*z[N])*integration(&Asypow::fN<N>,0,b);
    }
    ////////////////////////////////////////////////////////////////////////
    template<int L>
    double gL(double u){
        return z[L]*exp(theta*z[L])+a1(u)*exp(theta*z[L])+a2(u)*z[L]+a3(u);
    }
    template<int L>
    double fL(double t){
        return (b-t)*exp(-exp(theta*z[L])*lambda0*t)*pow( integration((&Asypow::gL<L>),0,t), 2 );
    }
    template<int L>
    double hL(double c){
        return exp(-exp(theta*z[L])*lambda0*c)*pow( integration(&Asypow::gL<L>,0,c), 2 );
    }

    template<int L>
    double termL()
    {
        const int N = 1000;
        double x1[N], x2[N], yf[N], yh1[N], yh2[N];
        for(int i=0;i<N;i++) {
            x1[i]=a*((double)(i+1))/((double)N);
            x2[i]=a+(b-a)*((double)(i+1))/((double)N);
            yf[i]=fL<L>(x2[i]);
            yh1[i]=hL<L>(x1[i]);
            yh2[i]=hL<L>(x2[i]);
        }
        double integralf = (b-a)*mean(yf,N);
        double integralh1 = a*mean(yh1,N);
        double integralh2 = (b-a)*mean(yh2,N);
        return p[L]*exp(theta*z[L])
              *(integralf +(b-a)*integralh1 +integralh2/(exp(theta*z[L])*lambda0));
    }
    ////////////////////////////////////////////////////////////////////
    // calculate the mean of eta
    template<int M>
    double gM(double u) {
        return (z[M]*exp(theta*z[M])+a1(u)*exp(theta*z[M])+a2(u)*z[M]+a3(u))
               *s_c(u)*exp(-exp(theta*z[M])*lambda0*u);
    }
    ////////////////////////////////////////////////////////////////////
    //covariance
    template<int K>
    double gK(double u) {
       return z[K]*exp(theta*z[K])+a1(u)*exp(theta*z[K])+a2(u)*z[K]+a3(u);
    }
    template<int K>
    double hK(double u) {
        return z[K]-e0(u);
    }
    template<int K>
    double fK(double t) {
        return (b-t)*(integration(&Asypow::gK<K>,0,t))
               *( z[K]-e0(t)-lambda0*exp(theta*z[K])*(integration(&Asypow::hK<K>,0,t)) )
               *exp(-exp(theta*z[K])*lambda0*t)
              -exp(-exp(theta*z[K])*lambda0*t)
               *(integration(&Asypow::gK<K>,0,t))
               *(integration(&Asypow::hK<K>,0,t));
    }
    template<int K>
    double mK(double t) {
           return (integration(&Asypow::gK<K>,0,t))*(z[K]-e0(t)-lambda0*exp(theta*z[K])
                 *(integration(&Asypow::hK<K>,0,t)))*exp(-exp(theta*z[K])*lambda0*t);
    }

    template<int K>
    double termK()
    {
        const int N = 1000;
        double d[N-1], x[N-1], yf[N-1], ym[N-1];
        for(int i=0;i<N-1;i++) {
            d[i]=a*((double)(i+1))/((double)N);
            x[i]=a+(b-a)*((double)(i+1))/((double)N);
            yf[i]=fK<K>(x[i]);
            ym[i]=mK<K>(d[i]);
        }
        double integralf = (b-a)*mean(yf,N-1);
        double integralm = a*mean(ym,N-1);
        return p[K]*exp(theta*z[K])*(integralf +(b-a)*integralm);
    }
    ///////////////////////////////////////
    Rcpp::NumericVector run()
    {
        double mu = sqrt(n)*lambda0*(integration(&Asypow::s1,0,b)-integration(&Asypow::f,0,b));
        double I0 = integration(&Asypow::g,0,b);
        double v1 = lambda0 *( vN<0>() + vN<1>() + vN<2>() );

        double power0= 1-pnormStd( ( sqrt(I0)*qnormStd(1-alpha/2)-mu)/sqrt(v1) )
                        +pnormStd( (-sqrt(I0)*qnormStd(1-alpha/2)-mu)/sqrt(v1) );


        double vapprox= v1;
        double v2;
        double power;
        double cov12;
        double vexact;
        double SSE;
        double variance;
        double v12;
        double diff;

        if(exactvar) {
            double Eeta = lambda0*p[0]*integration(&Asypow::gM<0>,0,b)
                         +lambda0*p[1]*integration(&Asypow::gM<1>,0,b)
                         +lambda0*p[2]*integration(&Asypow::gM<2>,0,b);
            v2    = (termL<0>() + termL<1>() + termL<2>())*lambda0*lambda0*lambda0/(b-a);
            v2 = v2-Eeta*Eeta; //variance of eta
            cov12 = (termK<0>() + termK<1>() + termK<2>())*lambda0*lambda0/(b-a);
            variance = v1+v2+2*cov12;
            power = 1-pnormStd( ( sqrt(I0)*qnormStd(1-alpha/2)-mu)/sqrt(variance) )
                                +pnormStd( (-sqrt(I0)*qnormStd(1-alpha/2)-mu)/sqrt(variance) );
            vexact = v1+v2+cov12;
            SSE = vapprox/vexact; 
            v12  =cov12;
            diff = v2+cov12;
        }
        else {
            v2       = NA_REAL;
            cov12    = NA_REAL;
            variance = NA_REAL;
            power    = NA_REAL;
            vexact   = NA_REAL;
            SSE      = NA_REAL;
            v12      = NA_REAL;
            diff     = NA_REAL;
        }

        return Rcpp::NumericVector::create( Rcpp::_["power"  ] = power,
                                            Rcpp::_["power0" ] = power0,
                                            Rcpp::_["v1"     ] = v1,
                                            Rcpp::_["v2"     ] = v2,
                                            Rcpp::_["v12"    ] = v12,
                                            Rcpp::_["vapprox"] = vapprox,
                                            Rcpp::_["vexact" ] = vexact,
                                            Rcpp::_["diff"   ] = diff,
                                            Rcpp::_["SSE"    ] = SSE );
    }
    ////////////////////////////////////
    double integration(double (Asypow::*mFunc)(double), double a, double b)
    {
        FunctionMaker<Asypow>::a= this;
        FunctionMaker<Asypow>::func=mFunc;

        return gslIntegration(FunctionMaker<Asypow>::funcForIntegration,a,b);
    }
};

RcppExport SEXP asypowRcpp(SEXP _n, SEXP _theta, SEXP _a, SEXP _b, SEXP _lambda0, SEXP _q, 
                           SEXP _p, SEXP _alpha, SEXP _z, SEXP _exactvar) {
BEGIN_RCPP
    Asypow asypow(Rcpp::as<double>( _n),
                  Rcpp::as<double>( _theta),
                  Rcpp::as<double>( _a),
                  Rcpp::as<double>( _b),
                  Rcpp::as<double>( _lambda0),
                  Rcpp::as<double>( _q),
                  Rcpp::as<Rcpp::NumericVector>(_p),
                  Rcpp::as<double>( _alpha),
                  Rcpp::as<Rcpp::NumericVector>(_z),
                  Rcpp::as<bool>(_exactvar) );
    return asypow.run();
END_RCPP
}


