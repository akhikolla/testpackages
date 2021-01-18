#include "../inst/include/entropyExport.h"


#include<Rcpp.h>
using namespace Rcpp;


// 'The entropy functions
RCPP_MODULE (moduleInfo) {
  function( "entropy_disc", &entropy_disc, "discrete entropy");
  function( "mutualInformation_disc", &mutualInformation_disc, "discrete MI");
  function( "mutualInformation_disc_u", &mutualInformation_disc_u, "discrete MI");
  function( "transferEntropy_disc", &transferEntropy_disc, "discrete TE");
  function( "entropy_cont", &entropy_cont, "Continuous entropy");
  function( "mutualInformation_cont", &mutualInformation_cont, "Cont MI");
  function( "transferEntropy_cont", &transferEntropy_cont, "Cont TE");
}
