#ifndef CLAMSEG_H_SMALLSCALES
#define CLAMSEG_H_SMALLSCALES

#include <list>

#include <Rcpp.h>

using namespace Rcpp;

/***************
* class small scales
* maintains additional finding on small scales
* Florian Pein, 2017
***************/
class SmallScales {
  private:
    unsigned int left_;
    unsigned int right_;
    unsigned int li_;
    unsigned int ri_;
    double stat_;
    bool noDeconvolution_;
  
  public:
    static std::list<SmallScales> listSmallScales_;
    static std::list<SmallScales>::iterator it_;
    
    SmallScales();
    SmallScales(unsigned int left, unsigned int right, unsigned int li, unsigned int ri,
                double stat, bool noDeconvolution);
    
    unsigned int left();
    unsigned int right();
    unsigned int li();
    unsigned int ri();
    double stat();
    bool noDeconvolution();

    static void update(unsigned int start, unsigned int len, double stat);
    void replace(unsigned int start, unsigned int len, unsigned int li, unsigned int ri,
                 double stat, bool noDe);
    void extend(unsigned int li, unsigned int ri);
    
    // clean up of static variables
    static void cleanUpGlobalVariables();
};

#endif
