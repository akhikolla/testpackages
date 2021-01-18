#include "ConConPiWiFun.h"


//RCPP_MODULE(mod_cplfunctionR){
//  using namespace Rcpp;
//
//	class_<cplfunctionR>( "cplfunctionR" )
//	//constructors
//	.constructor()
//	.constructor<double,double>()
//	//.constructor<double,double,double>()
//	.constructor<double,double,double,double>()
//	.constructor<Rcpp::NumericVector,Rcpp::NumericVector,double>()
//	;
//}

RCPP_MODULE(mod_cplfunction){
  using namespace Rcpp;

	class_<cplfunction>( "cplfunction" )
	//constructors
	.constructor()
	.constructor<double,double>()
	//.constructor<double,double,double>()
	.constructor<double,double,double,double>()
	.constructor<Rcpp::NumericVector,Rcpp::NumericVector,double>()

	.method("clone", &cplfunction::clone)
	//.field_readonly( "Breakpoints_", &cplfunction::get_BreakPoints_ )
	.field( "FirstBreakVal_", &cplfunction::FirstBreakVal_ )

	//methods
	.method("get_BreakPoints_",&cplfunction::get_BreakPoints)
	.method("Argmin",&cplfunction::Argmin)
	.method("Squeeze",&cplfunction::Squeeze)
	.method("Swap",&cplfunction::Swap)
	.method("Etoile",&cplfunction::Etoile)
	.method("eq",&cplfunction::eq)
	.method("Legendre",&cplfunction::Legendre)
	.method("EpiSum_Withline",&cplfunction::EpiSum_Withline)
	//.method("flip_push_left_pull_right",&cplfunction::flip_push_left_pull_right)
	//.method("flip_push_left",&cplfunction::flip_push_left)
 // .finalizer( &finalizer_of_cplfunction)
	;

  class_<cplfunctionvec>( "cplfunctionvec")
  .constructor()
  .constructor<int>()
  .method( "size", &cplfunctionvec::size)
 // .method("capacity", &cplfunctionvec::capacity,"Return size of allocated storage capacity. Returns the size of the storage space currently allocated for the vector, expressed in terms of elements.")
//  .method( "max_size", &cplfunctionvec::max_size)
  .method( "push_back", &cplfunctionvec::push_back )
//  .const_method( "at", &cplfunctionvec::at )
  .method("[[",&cplfunctionvec::vec_get)
  .method("[[<-",&cplfunctionvec::vec_set)
  .method("OptimMargInt",&cplfunctionvec::OptimMargInt,"Solves optimisation problem")
 // .method("OptimMargInt2",&cplfunctionvec::OptimMargInt2,"Solves optimisation problem")
  .method("SerialPush_1Breaks_Functions",&cplfunctionvec::SerialPush_1Breaks_Functions)
  .method("SerialPush_2Breaks_Functions",&cplfunctionvec::SerialPush_2Breaks_Functions)
  .method("SerialPush_nBreaks_Functions",&cplfunctionvec::SerialPush_nBreaks_Functions)
  .method("SerialPenalize",&cplfunctionvec::SerialPenalize)
  .method("SerialPush_Store_Functions",&cplfunctionvec::SerialPush_Store_Functions)
  ;
  
    Rcpp::function("OptimPriceMarket_l",&OptimPriceMarket_l)
    ;
  Rcpp::function("OptimPriceStorage_",&OptimPriceStorage_)
  ;


  Rcpp::function("SerialOptimPriceStorage",&SerialOptimPriceStorage)
  ;
  Rcpp::function("Suml",&Suml,"This function allows to sum two functions of class Rcpp_cplfunction. It does not modify the imput functions.")
  ;
  Rcpp::function("InfConvl",&InfConv,"This function performs infimum convolution of two functions of class Rcpp_cplfunction.")
  ;
  Rcpp::function("OptimPriceMarket_l",&OptimPriceMarket_l,"This function compute the market merit order.")
    ;

}


RCPP_MODULE(mod_cpqfunction){
  using namespace Rcpp;

	class_<cpqfunction>( "cpqfunction" )
	//constructors
	.constructor()
//	.constructor<double,double>()
	//.constructor<double,double,double>()
	//.constructor<double,double,double,double>()
	.constructor<Rcpp::NumericVector,Rcpp::NumericVector,Rcpp::NumericVector,double>()

	.method("clone", &cpqfunction::clone)
	//.field_readonly( "Breakpoints_", &cpqfunction::get_BreakPoints_ )
	.field( "FirstBreakVal_", &cpqfunction::FirstBreakVal_ )

	//methods
	.method("get_BreakPoints_",&cpqfunction::get_BreakPoints)
	.method("Argmin",&cpqfunction::Argmin)
	.method("Squeeze",&cpqfunction::Squeeze)
	.method("Swap",&cpqfunction::Swap)
	.method("Etoile",&cpqfunction::Etoile)
	.method("eq",&cpqfunction::eq)
	.method("evalf",&cpqfunction::evalf)

 // .finalizer( &finalizer_of_cplfunction)
	;

  Rcpp::function("Sumq",&Sumq,"This function allows to sum two functions of class Rcpp_cpqfunction. It does not modify the imput functions.")
  ;
  Rcpp::function("InfConvq",&InfConvq,"This function performs infimum convolution of two functions of class Rcpp_cplfunction.")
  ;

  class_<cpqfunctionvec>( "cpqfunctionvec")
  .constructor()
  .constructor<int>()
  .method( "size", &cpqfunctionvec::size)
 // .method("capacity", &cpqfunctionvec::capacity,"Return size of allocated storage capacity. Returns the size of the storage space currently allocated for the vector, expressed in terms of elements.")
//  .method( "max_size", &cpqfunctionvec::max_size)
  .method( "push_back", &cpqfunctionvec::push_back )
//  .const_method( "at", &cpqfunctionvec::at )
  .method("[[",&cpqfunctionvec::vec_get)
  .method("[[<-",&cpqfunctionvec::vec_set)
  .method("OptimMargInt",&cpqfunctionvec::OptimMargInt,"Solves optimisation problem")
  .method("SerialPush_1Breaks_Functions",&cpqfunctionvec::SerialPush_1Breaks_Functions)
  .method("SerialPush_0Breaks_Functions",&cpqfunctionvec::SerialPush_0Breaks_Functions)
  //.method("OptimPriceMarket",&cpqfunctionvec::OptimPriceMarket)
  ;

  Rcpp::function("OptimPriceMarket_q",&OptimPriceMarket_q)
  ;




}
