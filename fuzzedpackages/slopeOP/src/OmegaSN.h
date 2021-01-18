// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef OMEGASN_H
#define OMEGASN_H

#include <math.h>
#include <vector>
#include <list>

#ifndef DISABLE_RCPP
#include "Rcpp.h"
#endif//DISABLE_RCPP

#include "math.h"

class OmegaSN
{
public:
  OmegaSN(std::vector< double >& values, double firstdata, unsigned int nbSeg, unsigned int n);
  ~OmegaSN();

  std::vector< int > GetChangepoints() const;
  std::vector< double > GetParameters() const;
  double GetGlobalCost() const;
  double GetPruning() const;

  double** preprocessing(std::vector< double >& data) const;
  void Q0init(std::vector< double >& data) const;

  void algoNULL(std::vector< double >& data);
  void algoISOTONIC(std::vector< double >& data);

  void backtracking(unsigned int n);

  ///////////////

private:
  unsigned int nbSegments;
  unsigned int nbStates;
  double* states;

  double** S12P; //vectors S1, S2 et SP
  double*** Q;
  unsigned int*** lastIndState;
  unsigned int*** lastChpt;

  std::vector< int > changepoints; ///vector of changepoints build by fpop (first index of each segment). size c
  std::vector< double > parameters; ///vector of means build by fpop. size c
  double globalCost;
  double pruning; /// between 0 and 1. 1 if no pruning.
};

#endif // OMEGASN_H

