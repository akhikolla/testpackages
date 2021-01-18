#include <fstream> // problem under Mac (!) : fstream must be included before R.h...
#include <numeric> //for inner_product
#include <limits>
#include <ctime>
#include <algorithm> // stable_sort
#include "Krigtypes.h"
#define R_NO_REMAP
#include "smooth.h" // (includes qr.h which includes R.h)
#include "smoothFriends.h"  // // FRIENDs de la class CSmooth

CSmooth* test;
CQR<CQRtype> QR_zut;
int fnevalcounter;
// one dimensional equation solver (bisection_search) and minimiser (brent)

using namespace std;

int safeprint(std::string somestring) {
#ifdef NO_R_CONSOLE
   std::cout<<somestring;
#else
   REprintf(somestring.c_str());
#endif
return(0);
}

int safeprint(long double somelongdouble) {
#ifdef NO_R_CONSOLE
   std::cout<<somelongdouble;
#else
   std::stringstream stst;
   stst<<double(somelongdouble)<<" ";  // problems if not cast to souble
   REprintf(stst.str().c_str());
#endif
return(0);
}

template <typename someType>
int safeprint(someType somenum) {
#ifdef NO_R_CONSOLE
   std::cout<<somenum;
#else
   std::stringstream stst;
   stst<<somenum<<" ";
   REprintf(stst.str().c_str());
#endif
return(0);
}


covTypedef bisection_search(covTypedef (*func)(covTypedef),covTypedef x1,covTypedef x2) {
// Pour que (*func) puisse acceder des variables membres d'une classe
// il faut d['e]clarer (*func) comme une fonction friend
// et dans la-dite fonction il faut pouvoir affecter les membres ? une instantiation donn['e]e (eg instantiationPtr->membre)
// pour cela j'utilise un pointeur instantiationPtr vers cette instantiation, d['e]fini au niveau global
// ce ptr doit ?tre d['e]clar['e] apr[`e]s la lecture de Krig.h, il est donc d['e]clar['e] dans Krig.cpp
// et extern dans Krig.h
// on ne peut le mettre dans Krigmain.cpp/Krigmain.h car Krigmain.cpp -> Krig.h ->qr.h ->Krigmain.h
// et alors l'extern dans Krigmain.h serait lu avant la lecture de d['e]claration de classe dans Krig.h
// (NB qr.h pas tr[`e]s d['e]placable car la d['e]claration de classe dans Krig.h fait r['e]f ? la classe CQR).
// Et voila le pourquoi du ptr and Krig.cpp et des fonctions friend.
/** Rq a posteriori: headers probab mal organis'es, plus sous utilisation de extern?
     par contre l'utilisation de friend semble OK. Cf bootstrap dans Genepop pour exemple plus simple. **/
//
//R: function (x1, x2, f, tol = 1e-07, niter = 25, f.extra = NA, upcross.level = 0)
// rtbis, Numerical recipes
//Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be re...ned until its accuracy is ...xacc.
    const covTypedef xacc=numeric_limits<covTypedef>::epsilon()*(abs(x1)+abs(x2))/2; //suggestion Num. Rec.
    const int maxit=-2*int(log(numeric_limits<covTypedef>::epsilon())/log(2.));
    int j;
    covTypedef dx,f,fmid,xmid,rtb;
    f=covTypedef((*func)(x1));
    fmid=covTypedef((*func)(x2));
    if (f*fmid >= 0.0) {
#ifdef NO_R_CONSOLE
    std::cerr<<"(!) From CSmooth::bisection_search() : Root must be bracketed for bisection. ";
    if (batchDebug) cin.get();
     exit(-1);
#else
     REprintf("(!) From CSmooth::bisection_search() : Root must be bracketed for bisection. \n");
#endif
        }
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so that f>0
/*std::cout<<x1<<" ";
std::cout<<x2<<" ";
std::cout<<dx<<" ";
std::cout<<rtb;getchar();*/
    for (j=1;j<=maxit;j++) { //lies at x+dx.
    fmid=covTypedef((*func)(xmid=rtb+(dx *= 0.5))); //Bisection loop.
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    if (f*fmid >= 0.0) {
#ifdef NO_R_CONSOLE
           std::cerr<<"(!) From CSmooth::bisection_search() : Too many bisections. ";
           if (batchDebug) cin.get();
           exit(-1);
#else
           Rf_error("(!) From CSmooth::bisection_search() : Too many bisections. \n");
#endif
        }
return(numeric_limits<covTypedef>::signaling_NaN());
}



//==============================================================================
// MEMBRES

CSmooth::CSmooth(Cpointls & ptls,double GCV, int verbosity) {
  allocated=false; //imp!
  QR_Tallocated=false; //imp!
  this->verbosity=verbosity;
  lambdaEst=std::numeric_limits<covTypedef>::quiet_NaN(); //imp!
  CovFnParam.resize(0);
  initCovFnParam.resize(0);
  if(GCV==0) {
    fgcvPtr = &Krig_fgcv<covTypedef>;
  } else {
    fgcvPtr = &match_fs2hat_pure_error<covTypedef>;
  }

  covFamily="Matern";
#ifdef NO_R_CONSOLE
  std::cout<<"Covariance family = Matern"<<std::endl;
#else
  if (verbosity>1) REprintf("Covariance family = Matern\n");
#endif
  xyFileName=ptls.getname();
  xy=ptls.getxy();
  sort_compress();
  if (!allocated) {
    axialArray=new covTypedef**[KgPtNbr];
    for ( int i=0; i<KgPtNbr; i++) {
      axialArray[i]=new covTypedef*[KgPtNbr];
      for ( int j=0; j<KgPtNbr; j++) {
        axialArray[i][j]=new covTypedef[fittedparamnbr];
      }
    }
    euclArray=new covTypedef*[KgPtNbr];
    for ( int i=0; i<KgPtNbr; i++) {euclArray[i]=new covTypedef[KgPtNbr];}
    //            invA_y=new covTypedef[KgPtNbr];
    axialFocal=new covTypedef*[KgPtNbr];
    for ( int i=0; i<KgPtNbr; i++) {axialFocal[i]=new covTypedef[fittedparamnbr];}
    euclFocal=new covTypedef[KgPtNbr];
    allocated=true;
  }
  if (allocated) {
    long double tmp;
    for ( int i=0; i<KgPtNbr; i++) {
      for ( int j=0; j<i; j++) {
        for (int k=0;k<fittedparamnbr;k++) {tmp=xy[i][k]-xy[j][k];axialArray[i][j][k]=tmp*tmp;}
      }
    }
  }
  yobs.resize(KgPtNbr);
  for (int ii=0;ii<KgPtNbr;ii++) {
    yobs[ii]=xy[ii][fittedparamnbr];
  }
  nFixef=1; // ordinary kriging by default
  hglmOffset=0;
}



int CSmooth::sort_compress() {
//   long double SS;
   KgPtNbr_with_repl=xy.size();
   if (int(xy[0].size())!=(fittedparamnbr+1)) {
#ifdef NO_R_CONSOLE
  		std::cerr << "(!) CSmooth::sort_compress() called on data with suspicious number of columns" << std::endl;
  		std::cerr << "xy[0].size() = "<< xy[0].size() << " fittedparamnbr+1 = "<<fittedparamnbr+1 << std::endl;
  		if (batchDebug) std::cin.get();
  		exit(-1);
#else
  		Rf_error("(!) CSmooth::sort_compress() called on data with suspicious number of columns");
#endif
   }
//std::cout<<xy[0][0]<<" "<<xy[0][1]<<" "<<xy[0][2];getchar();
   stable_sort(xy.begin(),xy.end(),compareX<ioType>); //does what we need for the next compression step
//std::cout<<std::endl<<"! "<<xy[0][0]<<" "<<xy[0][1]<<" "<<xy[0][2];getchar();

   hglm_y.resize(0);
   std::vector<std::vector<ioType> >::iterator xyit;
   for(xyit=xy.begin();xyit!=xy.end();xyit++) hglm_y.push_back(xyit->back()); // keeps orginal (uncompressed) y values for SORTED data
   /***  will now construct the compressed data ***/
   std::vector<std::vector<ioType> >uniquex_y(xy);
   W.resize(0);
   W2.resize(0);
//   SSv.resize(0);
   for (std::vector<std::vector<ioType> >::iterator it=uniquex_y.begin();it!=uniquex_y.end();it++) {
       it->pop_back(); //removes y column
   }
   int count;
   ioType ysum;
   long double y2sum,tmp;
   ioType ysumsum=0; // hglm only
   long double y2sumsum=0; // hglm only
   int uniqueCol=0; // hglm only
   pureSS=0;

//NB. Dans fields, pure.ss est bien ce que je calcule: la sum_uniques sum_{y ; x=unique}(y-my(y))^2
//MSE est SS divis['e] par nombre total de points en excedent de 'uniques'
//shat.pure.error est sqrt(MSE)
   xyit=xy.begin();
   // for hglm computations
   //KgPtNbr is not yet known, so we cannot push_back vectors of size KgPtNbr instead we create
   vector<int> ijvector(0); // ith value contains corresponding j of the nozero element of the Zmatrix
   for (std::vector<std::vector<ioType> >::iterator it=uniquex_y.begin();it!=(uniquex_y.end());it++) {
       count=1;
       // these sums are only for a single unique x ! :
       ysum=xyit->back(); // back = y value, de xyit pointant sur la premiere ligne d'une serie de coord potentiellement repetees
       y2sum=ysum*ysum;
       ijvector.push_back(uniqueCol); //always at least one such line...
       // now add stuff for each additional line of redundant x
       while ((it+count)!=uniquex_y.end() && (*it)==(*(it+count))) {//['e]galit['e] des x vectors)
          tmp=(xyit+count)->back();
          ysum+=tmp;
          y2sum+=tmp*tmp;
          count++;
          ijvector.push_back(uniqueCol);
       }
       uniqueCol++;
/*       if (count==1) {
//          SSv.push_back(0.);
       } else {
          tmp=y2sum-ysum*ysum/count; //SS, not MS... ysum*ysum/count=count*(ymean)^2
//          SSv.push_back(tmp);
          pureSS+=tmp;
       }*/
       if (count>1) pureSS+=y2sum-ysum*ysum/count; //SS, not MS... ysum*ysum/count=count*(ymean)^2
       // RHS is an estimate both of parametric sigma^2 and of (count=2) times the residual error
       ysumsum+=ysum;
       y2sumsum+=y2sum;
       W.push_back(count);
       //Rprintf("ici %u",count);
       W2.push_back(sqrt(covTypedef(count)));
       it->push_back(ysum/count);
       if (count>1) uniquex_y.erase(it+1,it+count);
       //then it++ will point to the next different x
       //but we need to increment the xyit pointer
       xyit+=count;
   }
//   ostream_vec_vector(uniquex_y,std::cout);
//   ostream_vector(xyweights,std::cout);
   xy=uniquex_y;
   KgPtNbr=xy.size(); // number of unique coordinates (new size after compression)
   Zmatrix.resize(KgPtNbr_with_repl);
   for (unsigned int ii=0;ii<KgPtNbr_with_repl;ii++) {
       Zmatrix[ii].resize(0);
       Zmatrix[ii].resize(KgPtNbr,0);
       Zmatrix[ii][ijvector[ii]]=1;
   }
   if (verbosity && int(KgPtNbr_with_repl)==KgPtNbr) {
#ifdef NO_R_CONSOLE
  		std::cerr << "CSmooth::sort_compress() detected no replicate y values for any X values" << std::endl;
  		std::cerr << "  Estimation of correlation parameters will fail. " << std::endl;
#else
      // Rf_error("(!) CSmooth::sort_compress() detected no replicate y values for any X values. I exit");
       REprintf("CSmooth::sort_compress() detected no replicate y values for any X values.\n");
#endif
   }
   pure_error=pureSS/(KgPtNbr_with_repl-KgPtNbr);
   hglmPhi=pure_error;
   phiFix=-1; // default, means that phi must be estimated by HLCor
   // more initial values for hglm algo
   hglmLambda=((y2sumsum-ysumsum*ysumsum/KgPtNbr_with_repl)-pureSS)/KgPtNbr; // let's forget any unbiasing correction here...
   lamFix=-1; // default, means that lambda must be estimated by HLCor
//cout<<hglmLambda<<" "<<pure_error;
   hglmBeta=ysumsum/KgPtNbr_with_repl; // also suboptimal for uncorrelated data but...
   //laborious code
   KgLow.resize(fittedparamnbr,numeric_limits<ioType>::max());
   KgUp.resize(fittedparamnbr,-numeric_limits<ioType>::max()); //numeric_limits<double>::min() is not what you think...
   for (unsigned int ii=0;ii<xy.size();ii++) {
       for (int jj=0;jj<fittedparamnbr;jj++) {
           if (xy[ii][jj]<KgLow[jj]) KgLow[jj]=xy[ii][jj];
           if (xy[ii][jj]>KgUp[jj]) KgUp[jj]=xy[ii][jj];
       }
   }
#ifdef NO_R_CONSOLE
   if (verbosity) {
     std::cout<<"\nLower bounds of Kriged points:\n";
     for (std::vector<ioType>::iterator ii=KgLow.begin();ii<KgLow.end();ii++)
         std::cout<<(*ii)<<" ";
     std::cout<<"\nUpper bounds of Kriged points:\n";
     for (std::vector<ioType>::iterator ii=KgUp.begin();ii<KgUp.end();ii++)
         std::cout<<(*ii)<<" ";
     std::cout<<"\nPure RMSE estimate from smoothed duplicated points: "<<sqrt(pure_error);
   }
#else
   if (verbosity) {
     Rprintf("\nLower bounds of Kriged points:\n");
     std::stringstream stst;
     std::string st="";
     for (std::vector<ioType>::iterator ii=KgLow.begin();ii<KgLow.end();ii++) {
        stst<<(*ii);
        st+=stst.str()+" ";
        stst.str("");
     }
     Rprintf("%s\n",st.c_str());
     st.resize(0);
     Rprintf("Upper bounds of Kriged points:\n");
     for (std::vector<ioType>::iterator ii=KgUp.begin();ii<KgUp.end();ii++){
        stst<<(*ii);
        st+=stst.str()+" ";
        stst.str("");
     }
     Rprintf("%s\n",st.c_str());
     st="Pure RMSE estimate from kriged duplicate points: ";
     stst<<sqrt(pure_error);
     st+=stst.str()+"\n";
     stst.str("");
     Rprintf("%s\n",st.c_str());
   }
#endif
return 0;
}


int CSmooth::filleuclArray() { // this function call only in Krig_engine_default
//std::cout<<"CovTheta2 au debut filleuclArray ";
//ostream_vector(CovTheta2,std::cout);
  long double euclidian;
  for (int ii=0;ii<KgPtNbr;ii++) {
     euclArray[ii][ii]=0;
     for (int jj=0;jj<ii;jj++) {
        euclidian=0;                                     /*axialArray contains *squared* axial dist...*/
        for (int kk=0;kk<fittedparamnbr;kk++) euclidian+=(axialArray[ii][jj][kk]/CovTheta2[kk]);
        euclArray[ii][jj]=sqrt(euclidian);
        if (ISNAN(euclidian) || ! R_FINITE(euclidian)) {
#ifdef NO_R_CONSOLE
            std::cerr<<"(!) From CSmooth::filleuclArray(): Inf/NaN euclidian distance"<<std::endl;
            std::cerr<<"pair i,j: "<<ii<<","<<jj<<std::endl<<"CovTheta: ";
            ostream_vector(CovTheta,std::cerr);
            if (batchDebug) cin.get();
            exit(-1);
#else
            std::stringstream stst;
            std::string st="";
            stst << "Squared axial distances:";
            for (int kk=0;kk<fittedparamnbr;kk++) {
              stst << axialArray[ii][jj][kk];
              st+=stst.str()+" ";
              stst.str("");
            }
            stst << "; CovTheta2:";
            for (int kk=0;kk<fittedparamnbr;kk++) {
              stst << CovTheta2[kk];
              st+=stst.str()+" ";
              stst.str("");
            }
            REprintf("%s\n",st.c_str());
            Rf_error("(!) From CSmooth::filleuclArray(): Inf/NaN euclidian distance");
#endif
        }
     }
   }
return 0;
}

int CSmooth::fillaxialFocal(std::vector<ioType>& focal) { //& useful
    ioType tmp; // hmmmmmmmmmmmmmmm
    for (int j=0; j<KgPtNbr; j++) {
        for (int k=0;k<fittedparamnbr;k++) {
/*#ifdef NO_R_CONSOLE
//            std::cout<<focal[k]<<" "<<xy[j][k];getchar();
#else
            REprintf("%d %d %Lf %Lf\n",j,k,focal[k],xy[j][k]);
#endif*/
            tmp=focal[k]-xy[j][k];
            axialFocal[j][k]=tmp*tmp;
        }
    }
return 0;
}


int CSmooth::filleuclFocal() {
  long double euclidian;

  for (int jj=0;jj<KgPtNbr;jj++) {
        euclidian=0;
        for (int kk=0;kk<fittedparamnbr;kk++) {
            //axialFocal[][] is a squared distance hence /theta^2
            // theta being a scale factor for a coordinate; CovTheta2 is already squared
          euclidian+=(axialFocal[jj][kk]/CovTheta2[kk]);

#ifdef NO_R_CONSOLE
#else
//                    REprintf("%f %f %f\n",axialFocal[jj][kk],pow(CovTheta[kk],covTypedef(2.)));
#endif
        }
        euclFocal[jj]=sqrt(euclidian);
     }

return 0;
}








