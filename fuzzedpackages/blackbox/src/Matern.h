#ifndef H_MATERN
#define H_MATERN

#include <cstdlib>
#include <string>
#include <sstream>
#include <limits>
#include <cmath> //for std::abs
#include "Rmath.h" // bessel_k (and apparently required before loading Rcpp.h...)
#include "Bessel_nr.h" // bessk
#include <Rcpp.h> // a charger apres les headers contenant des templated functions ???
//                  (sinon par exemple il ne comprend plus error())


namespace NS_GG {
 extern int a;
 extern covTypedef b;
}

template <typename Typeforcov>
Typeforcov gammln(Typeforcov xx) {
  Typeforcov xxx,tmp,ser;
  static Typeforcov cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.53363382e-5};
  xxx=xx-1.0;
  tmp=xxx+5.5;
  tmp-=(xxx+0.5)*log(tmp);
  ser=1.0;
  for(int j=0;j<5;j++)
  {xxx+=1.0;
    ser+=cof[j]/xxx;
  }
  return -tmp+log(2.50662827465*ser);
}



template <typename Typeforcov>
Typeforcov Matern(Typeforcov dist,const Typeforcov& smoothness) {
    /* use fields' parameterization, which is not Stein's one, p. 31*/
    Typeforcov cov,ddist;
    const Typeforcov ONE=1.;  // this looks a bit silly but (1) does not harm; (2) pow overloading....
    const Typeforcov TWO=2.;
    if (dist<std::numeric_limits<Typeforcov>::epsilon())
        return 1.;
    else {
        cov=pow(TWO,smoothness-ONE)*exp(gammln(smoothness));
        cov=1./cov;
        /** std::abs is not the same as abs (except for abs in context using namespace std;) **/
        if (false && std::fabs(smoothness-4)<std::numeric_limits<Typeforcov>::epsilon()) { // special case smoothness=4
            cov*=bessk<Typeforcov>(4,dist); // uses NR code for this integer case... seems OK, moderate gain in speed.
          // problem is that this is not continuous 'enough' with the genera case, actually impeding optimisation => false added 2015/11
            ddist=dist*dist;
            cov*=ddist*ddist;
        } else {
            cov*=R::bessel_k(dist,smoothness,1.);
            cov*=pow(dist,smoothness);
        }
//DEBUG
//cov*=0.9999;
/*        if (ISNAN(cov) || ! R_FINITE(cov)) {
            std::stringstream stst;
            stst<<"Matern pb: "<<dist<<" "<<smoothness<<std::endl;
#ifdef NO_R_CONSOLE
            std::cout<<stst.str();
#else
            REprintf(stst.str().c_str());
#endif
        }*/
        return cov;
    }
}

#endif

