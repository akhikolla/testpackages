// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef OMEGAOP_H
#define OMEGAOP_H

#include <math.h>
#include <vector>
#include <list>

#ifndef DISABLE_RCPP
#include "Rcpp.h"
#endif//DISABLE_RCPP

#include "math.h"

class OmegaOP
{
  public:
    OmegaOP(std::vector< double >& values, double firstdata, double beta, unsigned int n);
    ~OmegaOP();

    std::vector< int > GetChangepoints() const;
    std::vector< double > GetParameters() const;
    double GetGlobalCost() const;
    double GetPruning() const;

    double** preprocessing(std::vector< double >& data) const;

    void algo(std::vector< double >& data);
    void algoChannel(std::vector< double >& data);
    void algoPruning(std::vector< double >& data);

    void backtracking(unsigned int n);

    ///////////////
    void algoISOTONIC(std::vector< double >& data);
    void algoUNIMODAL(std::vector< double >& data);
    void algoSMOOTHING(std::vector< double >& data, double minAngle);

  private:
    double penalty;
    unsigned int nbStates;
    double* states;

    double** S12P; //vectors S1, S2 et SP
    double** Q;
    unsigned int** lastIndState;
    unsigned int** lastChpt;

    std::vector< int > changepoints; ///vector of changepoints build by fpop (first index of each segment). size c
    std::vector< double > parameters; ///vector of means build by fpop. size c
    double globalCost;
    double pruning; /// between 0 and 1. 1 if no pruning.
};

#endif // OMEGAOP_H
