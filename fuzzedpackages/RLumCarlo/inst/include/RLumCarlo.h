// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Title:  RLumCarlo Header
// Author:  Sebastian Kreutzer, Geography & Earth Sciences, Aberystwyth University (United Kingdom),
// Contact: sebastian.kreutzer@aber.ac.uk
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef RLUMCARLO_H
#define RLUMCARLO_H
#include <RcppArmadillo.h>
  //announce calc_detat()
  double calc_deltat(arma::vec times);

  //set constants
  //set Boltzmann's constant
  static const double k_B = 8.617333262145e-05;
#endif


