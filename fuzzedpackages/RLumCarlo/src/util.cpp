// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Title:   util.cpp
// Purpose: Internal C++ helper functions
// Author:  Sebastian Kreutzer, Geography & Earth Sciences, Aberystwyth University (United Kingdom)
// Contact: sebastian.kreutzer@aber.ac.uk
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "RLumCarlo.h"

// calculate delta t including error treatment
double calc_deltat(arma::vec times){
  //make sure we crash fast
  //check length
  if(times.size() < 2) {
    throw std::range_error("[RLumCarlo Internal Error] The length of times cannot be smaller than 2!");

  }

  //check if values are equidistant
  arma::vec v_test = arma::floor(arma::diff(times));
  if (!std::equal(v_test.begin() + 1, v_test.end(), v_test.begin())){
    throw std::range_error("[RLumCarlo Internal Error] Non-equidistant elements in times are not supported!");

  }

  //finally calculate delta t
  double delta_t = times[1] - times[0];

  //return
  return(delta_t);

}
