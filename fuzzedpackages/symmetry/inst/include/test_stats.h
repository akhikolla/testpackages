#ifndef symmetry_test_stats
#define symmetry_test_stats

using namespace Rcpp;

double MOI_Cpp(const NumericVector& X, double k);
double MOK_Cpp(const NumericVector& X, double k);
double K2_Cpp(const NumericVector& X);
double K2U_Cpp(const NumericVector& X);
double NAI_Cpp(const NumericVector& X, double k);
double NAK_Cpp(const NumericVector& X, double k);
double BHI_Cpp(const NumericVector& X);
double BHK_Cpp(const NumericVector& X);
double CM_Cpp(const NumericVector& X);
double MI_Cpp(const NumericVector& X);
double MGG_Cpp(const NumericVector& X);
double B1_Cpp(const NumericVector& X);
double NAC1_Cpp(const NumericVector& X, double a);
double NAC2_Cpp(const NumericVector& X, double a);
double BHC1_Cpp(const NumericVector& X, double a);
double BHC2_Cpp(const NumericVector& X, double a);
double HM_Cpp(const NumericVector& X, double a);
double FM_Cpp(const NumericVector& X);
double KS_Cpp(const NumericVector& X);
double SGN_Cpp(const NumericVector& X);
double WCX_Cpp(const NumericVector& X);
double RW_Cpp(const NumericVector& X);
double BH2_Cpp(const NumericVector& X);
#endif
