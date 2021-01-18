#ifndef H_SMOOTH
#define H_SMOOTH

//#define TEST_OCV

#include <vector>
#ifdef NO_R_CONSOLE
#include <iostream>
#endif
#include <cmath>
#include <limits>
#include <ctime>
#include <numeric> //for inner_product
#include "qr.h" // (includes R.h)
#include "pointls.h"
#include <Rcpp.h>
#include <RcppEigen.h> //Eigen::SelfAdjointEigenSolve

#define CGOLD 0.381966011250105151795413165634361882279690820194237137864551 //(3-sqrt(5))/2
//ZEPS is a small number that protects against trying to achieve fractional accuracy for a minimum that
//happens to be exactly zero.
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

extern int fnevalcounter;

double cylindercov(double euclidian);
covTypedef bisection_search(covTypedef (*func)(covTypedef),covTypedef x1,covTypedef x2);

int safeprint(std::string somestring);
int safeprint(long double somelongdouble);
template <typename someType>
int safeprint(someType somenum);

template <typename Typeforcov>
Typeforcov match_fs2hat_pure_error(Typeforcov lambda);

template <typename Typeforcov>
Typeforcov Krig_fgcv(Typeforcov lambda);



// this does the same as the default < in stable_sort but
//(1) ignores the last element of the vectors
// detects if matrix was unsorted when called from R
template <typename someType>
bool compareX(std::vector<someType> a, std::vector<someType> b) {
    typename std::vector<someType>::iterator ita=a.begin();
    typename std::vector<someType>::iterator itb=b.begin();
//getchar();
//ostream_vector(a,std::cout);
//ostream_vector(b,std::cout);
// ita est plus loin dans la matrice que itb... si c'est bien tri['e],
// *ita>*itb, le retour est false et rien n'est retri['e]
    while (ita<(--a.end())) { // ignores the last element of the vectors !!
        if ((*ita)>(*itb)) {
            return(false);
        } else {
#ifndef NO_R_CONSOLE // NOT NO_R = if DLL
// called from R code: the arrya should already be sorted so that the order of
// krig coefficients is the same in C and R code
            if ((*ita)!=(*itb)) Rf_error("(!) From compareX() in DLL : parameter points provided by R call not sorted. \n");
#endif
            ita++;itb++;
        }
    }
    return true;
}

template <typename TypeforBrent>  /// internal type
covTypedef brent(covTypedef (*f)(covTypedef),covTypedef ax, covTypedef bx, covTypedef cx, covTypedef* xmin) {
const TypeforBrent EPSILON=(std::numeric_limits<TypeforBrent>::epsilon());
const TypeforBrent TOL=(sqrt(EPSILON));
const TypeforBrent ZEPS=(100.*EPSILON);
//Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
//between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
//the minimum to a fractional precision of about TOL using Brent's method. The abscissa of
//the minimum is returned as xmin, and the minimum function value is returned as brent, the
//returned function value.
///Brent's book p.53:
//maximum # of required [fn evaluations] is (k+1)?-2
//where k=intsup(Log_2([initial interval width]/(2*machine_epsilon|x|+t)))
// here t seems to be the absolute tolerance
//|x| where x is the abcissa (note miniminisation over the interval)
//cf computation of tol1 and tol2
//Brent observes that in no case more than 3(k+1) iterations were needed.
// j'ai pris 4 k ci dessous.
    //             [   ?           .           +|                                     |+     .?                  ]
    int ITMAX=4*int(log(std::abs(ax-cx)/(2.*EPSILON*((std::abs(ax) < std::abs(cx) ? std::abs(ax) : std::abs(cx)))+ZEPS))/log(2.)+1-EPSILON);
    if (ITMAX<4) ITMAX=4; // protection from silly case where ax-cx is small relative to epsilon => negative log/log....
    int iter;
    TypeforBrent a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    TypeforBrent ebrent=0.0; //This will be the distance moved on the step before last.
    a=(ax < cx ? ax : cx); //a and b must be in ascending order, but input abscissas need not be.
    b=(ax > cx ? ax : cx);
    x=w=v=bx; //Initializations...
    fw=fv=fx=covTypedef((*f)(x));
    for (iter=1;iter<=ITMAX;iter++) { //Main program loop.
        xm=0.5*(a+b);
        tol2=2.0*(tol1=TOL*std::abs(x)+ZEPS);
        if (std::abs(x-xm) <= (tol2-0.5*(b-a))) { //Test for done here.
            *xmin=x;
            return fx;
        }
        if (std::abs(ebrent) > tol1) { //Construct a trial parabolic fit.
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=std::abs(q);
            etemp=ebrent;
            ebrent=d;
            // here g++ (4.6.3) -03 or -O2 says that ebrent may be used uninitialized...
            if (std::abs(p) >= std::abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
               //ebrent=(x >= xm ? a-x : b-x);
               //if (x>=xm) {ebrent=a-x;} else {ebrent=b-x;}
               d=CGOLD*(ebrent=(x >= xm ? a-x : b-x));
        //    The above conditions determine the acceptability of the parabolic fit. Here we
        //    take the golden section step into the larger of the two segments.
            } else {
                d=p/q; //Take the parabolic step.
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                d=SIGN(tol1,xm-x);
            }
        } else {
            d=CGOLD*(ebrent=(x >= xm ? a-x : b-x));
        }
        u=(std::abs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=covTypedef((*f)(u)); //This is the one function evaluation per iteration.
//DEBUG
//{ std::stringstream stst;stst<<"      "<<u<<" "<<fu<<std::endl;REprintf(stst.str().c_str());}
        if (fu <= fx) { //Now decide what to do with our function evaluation.
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u) //Housekeeping follows:
            SHFT(fv,fw,fx,fu)
        } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu <= fv || v == x || v == w) {
                v=u;
                fv=fu;
            }
        } //Done with housekeeping. Back for another iteration.
    }
//    *xmin=x; //Never get here.
#ifdef NO_R_CONSOLE
    std::cerr<<"(!) From CSmooth::brent(): Too many iterations ("<<ITMAX<<").";
    if (batchDebug) std::cin.get();
    exit(-1);
#else
     REprintf("%d iterations.",ITMAX); // typically goes to
     Rf_error("(!) From CSmooth::brent(): Too many iterations.\n");
#endif
return fx;
}


class CSmooth { //class optimized for kriging operations only
    std::string xyFileName;
    std::vector<std::vector<ioType> > xy;
    std::vector<ioType>hglm_y; // keeps the response in uncompressed form
    // the hglm ranefs on two scales
    std::vector<ioType>u_h;
    std::vector<ioType>v_h;
    std::vector<std::vector<int> > Zmatrix; // for hglm computations
    std::vector<std::vector<covTypedef> >  ZL;
    covTypedef*** axialArray; //unscaled squared axial distances (only j<i half is used)
    covTypedef** euclArray; //euclidian scaled distance (only j<i half is used)
    covTypedef* euclFocal;
    std::vector<std::vector<covTypedef> >  covMat; // covariance matrix (only j<i half is used? No!)
    covTypedef** axialFocal; //unscaled squared axial distances
    std::vector<covTypedef> covFocal; // covariance vector
    int ncolT; // nb cols Tmatrix
    CQR<CQRtype> *QR_T;    // implements the QR decomp of T //use ptr in order not to call any constructor now.
    bool QR_Tallocated; // to control deallocation QR_T ptr at destruction CSmooth instantiation
    std::vector<ioType>  yobs;
    /*** ideally I would template the CSmooth class with a Typeforcov. But too big a change, not urgent.
         Anyway I should first ensure that the code is OK with covTypedef without template***/
    std::vector<internalTypeEigen>  u; // une des sorties essentielles de Krig.engine.default // but not simply i/o with R !
    // these coefficients are used together with the covariances in predict, hence their covType
    std::vector<internalTypeEigen> D_invEigVals; //idem : collage de 0 et de 1/eigenvalues
    std::vector<std::vector<internalTypeEigen> >  eigVecs; // idem : collage de 0 et d'eigenvectors
#ifdef TEST_OCV
    std::vector<std::vector<internalTypeEigen> >  TeigVecs2; //same as eigVecs but transposed and squared
#endif
    std::vector<std::vector<internalTypeEigen> >  G_for_Krig_c_coef; // idem : ....
//    internalTypeEigen* invA_y; // stores inv(A).y in predictor c.inv(A).y
    std::vector<covTypedef> CovTheta;
    std::vector<covTypedef> CovTheta2; // squared values to save some computation time
    covTypedef df; //degrees of freedom in gcv computations
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    std::vector<covTypedef> grid_df; //? stocker pour graphique ?
    covTypedef smoothness;
    covTypedef minSmoothness; // added 03/2011, to control min Smoothness in call from R
    int KgPtNbr;
    unsigned int KgPtNbr_with_repl; //nb "points" (lignes) avant compression := sum(W)
    std::string covFamily;
    bool allocated;
    std::vector<std::vector<covTypedef> >Tmatrix; //*vecteur* de 1 de la longueur des ptls
    //      QR<covTypedef> QR_T; //will contain QR decomp of Tmatrix
    std::vector<std::vector<covTypedef> > QT; //will contain Q from QR decomp of Tmatrix
    std::vector<std::vector<covTypedef> > RT; //will contain R from QR decomp of Tmatrix
    std::vector<covTypedef> gcv_grid;
    std::vector<covTypedef> lambda_grid;
    covTypedef lambdaEst; //gcv estimator of lambda; // !! lambda_GCV <=> phi_HGLM/lambda_HGLM
    std::vector<internalTypeEigen> coefs_fixed;
    std::vector<internalTypeEigen> coefs_random;
    std::vector<covTypedef> initCovFnParam; //20151024 non-default initial values of cov fn params (theta+smoothness)
    std::vector<covTypedef> CovFnParam; //uniquement pour stocker meilleures estimations
    covTypedef GCVmini;
    std::vector<int> W;
    std::vector<covTypedef> W2;
//    std::vector<long double> SSv; // vector of \sum y_i^2 -\sum \bar{y}^2 per unique coordinates
    long double pureSS;
    ioType hglmLambda;
    ioType hglmPhi;
    ioType lamFix;
    ioType phiFix;
    ioType hglmBeta;
    ioType hglmOffset;
    int verbosity;
    double p_bv;
    double logabsdetAugdesign;
    int nFixef;
    std::vector<ioType> hglmTheta;
    covTypedef pure_error; //MSE
    std::vector<ioType> KgLow;
    std::vector<ioType> KgUp;
    std::vector<ioType> maxrange;
    std::vector<bool> estimBools; // infos sur les cov params a estimer
    std::vector<covTypedef> lowerbound;
    std::vector<covTypedef> upperbound;
    covTypedef (*fgcvPtr)(covTypedef); //types must be the same as those given in the assignment of fgcvPtr
public:
      CSmooth() {allocated=false;QR_Tallocated=false;CovFnParam.resize(0);
               lambdaEst=std::numeric_limits<covTypedef>::quiet_NaN();
               xyFileName="noptls_0";}; //
      //The following constructor creates required stuff up to axial squared distances
      //Computation of covariances for the 'x' table best done through predictor()
      //Computation of covariances for focal points best done through predict()
       CSmooth(Cpointls & ptls,double GCV, int verbosity);
      ~CSmooth() { //new dans CSmooth(ioType** xy,...);
        if (allocated) {
            for (int i=0; i<KgPtNbr; i++) {
                for (int j=0; j<KgPtNbr; j++) {
                    delete[] axialArray[i][j];
                }
                delete[] axialArray[i];
            }
            delete[] axialArray;
            for (int i=0; i<KgPtNbr; i++) delete[] euclArray[i];
            delete[] euclArray;
//            delete[] invA_y;
            for (int i=0; i<KgPtNbr; i++) {
                delete[] axialFocal[i];
            }
            delete[] axialFocal;
            delete[] euclFocal;
            allocated=false;
        }
        if (QR_Tallocated) {delete QR_T;QR_Tallocated=false;}
      };
      int filleuclArray();

      template<typename Typeforcov>
      int fillcovMat(const Typeforcov& smoothness);

      int fillaxialFocal(std::vector<ioType>& focal);
      int filleuclFocal();

      template<typename Type>
      int fillcovFocal();

      template <typename Typeforcov>
      ioType predict(std::vector<ioType> focal,std::string method="QR");
//      int evaluatefixed();

      template <typename Typeforcov,typename Type>  // remarquablement plus lent si on passe CovTheta par reference !!
      int Krig_engine_default(std::vector<Typeforcov> CovTheta, const Typeforcov& smoothness);

      template <typename Typeforcov,typename Type>  // remarquablement plus lent si on passe CovTheta par reference !!
      int HLCor(std::vector<Typeforcov> covtheta, const Typeforcov& smoothness);

      template <typename Typeforcov,typename Type>  // remarquablement plus lent si on passe CovTheta par reference !!
      int augLinMod(std::vector<Typeforcov> covtheta, const Typeforcov& smoothness);

      template <typename Typeforcov>
      Typeforcov gcv_Krig(); //estimates lambda AND returns GCV value

      template <typename Typeforcov>
      std::vector<Typeforcov> logLik_hglm(std::vector<std::vector<Typeforcov> >* newZL); //estimates lambda AND returns GCV value

      template <typename Typeforcov>
      Typeforcov ocv_Krig(); //estimates lambda AND returns OCV value
//***********
      template <typename Typeforcov>
      friend Typeforcov Krig_fgcv(Typeforcov lambda);

      template <typename Typeforcov>
      friend Typeforcov match_fs2hat_pure_error(Typeforcov lambda);
#ifdef TEST_OCV
      template <typename Typeforcov>
      friend Typeforcov Krig_ocv(Typeforcov lambda);
#endif
//************
      template<typename Typeforcov>
      friend Typeforcov Krig_fdf(Typeforcov loglambda);

      template<typename Typeforcov>
      Typeforcov Krig_df_to_lambda(Typeforcov explicitdf=-1);

      template<typename Typeforcov,typename TypeforES>
      int Krig_coef(Typeforcov lambda=std::numeric_limits<Typeforcov>::quiet_NaN());

      template<typename Typeforcov,typename TypeforES>
      int hglm_Krig_coef();

      template <typename Typeforcov,typename TypeforES>
      Typeforcov GCV_lamVar_covFix(std::vector<Typeforcov> covparam);

      template <typename Typeforcov,typename TypeforES>
      friend Typeforcov optimGivenOffset(ioType locOffset);

      template <typename Typeforcov,typename TypeforES>
      void optimOverOffset(bool resetPredictor,double maxSmoothness,bool verbose=true);
      friend SEXP GCV_lamVar_covFix_Wrapper(SEXP a, SEXP fixedSmoothness, SEXP returnFnvalue);
      friend SEXP CcovFocal(SEXP focal,SEXP CKrigidxP);
      friend Rcpp::List Krig_coef_Wrapper(SEXP aA, SEXP lambdaP);
      friend void testfn(); //to access private member in test fn
      friend bool intern_newCSmooth( //friend to CSmooth declaration
          double *xy,
          int *nrowxy, //with replicates
          int *ncolxy,
          int *nuniquerows, // required to tell the allocated size of c, d, D, u in *R*
          //double *maxSmoothness,
          //    double *c,
          //    double *d,
          //    double *D,
          //    double *u,
          //    double *lambda,
          double *GCV,
          //    double *covfnparam, //p+1 form
          // int *fneval,
          int *optimiseBool,
          int *verbosity
          //    double *initcovfnparam //p+1 form  // new arg 2015/10/24
      );
      friend int deleteCSmooth();
      template <typename Typeforcov,typename TypeforES>
      void hglmjoint(bool resetPredictor, double maxSmoothness,bool verbose=true);
      template <typename Typeforcov,typename TypeforES>
      void optimOverCorrPars(bool resetPredictor, double maxSmoothness,bool verbose=true);

      int sort_compress();

};

extern CSmooth* test;

#endif
