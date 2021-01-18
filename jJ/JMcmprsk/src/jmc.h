#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>


#include <Rcpp.h>
using namespace Rcpp;

namespace jmcspace {

//const std::vector<double> xs;
//const std::vector<double> ws;


 double HAZ(const gsl_matrix *H, const double t);

 double CH(const gsl_matrix *H, const double t);

 double MulVV(const gsl_vector *Z,const gsl_vector *beta);

 void MulV(const gsl_vector *Z,gsl_matrix *ZZ);

 void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta);

 void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB);

 int inv_matrix(gsl_matrix *x_square);

 double Abs(const double a, const double b);

 int DiffM1(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

 double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

 double DiffV(const gsl_vector *veca, const gsl_vector *vecb);

 double Min(const double t1, const double t2);

 void STAT(gsl_matrix *store,int i,double *mean,double *sd);

 int DiffM2(const gsl_matrix *preH1,const gsl_matrix *H1,const gsl_matrix *preH2,const gsl_matrix *H2);

 void TransM(const gsl_matrix *A, gsl_matrix *B);

 int Sbeta(gsl_vector *beta, double *sigma, const gsl_matrix *Y, const int p1a);




 int EM(
       gsl_vector *beta,
       gsl_matrix *gamma,
       gsl_vector *vee,
       gsl_matrix *H01,
       gsl_matrix *H02,
       double *sigma,
       gsl_matrix *sig,
       const gsl_matrix *Y,
       const gsl_matrix *C,
       const gsl_vector *M1,
       const int p1a,
       const int maxl,
       const int point,
       const std::vector<double> xs,  
       const std::vector<double> ws
       );


 int GetCov(
           gsl_matrix *Cov,
           gsl_vector *beta,
           const gsl_matrix *gamma,
           const gsl_vector *vee,
           const gsl_matrix *H01,
           const gsl_matrix *H02,
           const double sigma,
           gsl_matrix *sig,
           const gsl_matrix *Y,
           const gsl_matrix *C,
           const gsl_vector *M1,
           const int p1a,
           const int maxl,
           const int point,
           const std::vector<double> xs,  
           const std::vector<double> ws
           );


 int GetE(
          gsl_vector *FUNU,
          gsl_vector *FUNUS,
          gsl_matrix *FUNB,
          gsl_matrix *FUNBS,
          gsl_matrix *FUNBU,
          gsl_matrix *FUNE,
          gsl_matrix *FUNUSE,
          gsl_matrix *FUNUE,
          gsl_matrix *FUNW,
          gsl_matrix *FUNWB,
          gsl_matrix *FUNWBS,
          const gsl_vector *beta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const double sigma,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int maxl,
          const int point,
          const std::vector<double> xs,  
          const std::vector<double> ws
          );



 double Getloglik(
          const gsl_vector *beta,
          const gsl_matrix *gamma,
          const gsl_vector *vee,
          const gsl_matrix *H01,
          const gsl_matrix *H02,
          const double sigma,
          const gsl_matrix *sig,
          const gsl_matrix *Y,
          const gsl_matrix *C,
          const gsl_vector *M1,
          const int p1a,
          const int maxl,
          const int point,
          const std::vector<double> xs,  
          const std::vector<double> ws
          );


 int Diff(
         const gsl_vector *prebeta,
         const gsl_vector *beta,
         const gsl_matrix *pregamma,
         const gsl_matrix *gamma,
         const gsl_vector *prevee,
         const gsl_vector *vee,
         const gsl_matrix *preH01,
         const gsl_matrix *H01,
         const gsl_matrix *preH02,
         const gsl_matrix *H02,
         const double presigma,
         const double sigma,
         const gsl_matrix *presig,
         const gsl_matrix *sig
         );
//declare before use it
Rcpp::List jmc_cmain(int k, int n1,int p1,int p2, int maxl, int p1a, int maxiterations, int point, std::vector<double>xs,  std::vector<double> ws, std::string yfile, std::string cfile, std::string mfile, int  trace);		 
#define array_size 100000
#define kappa 100
}
