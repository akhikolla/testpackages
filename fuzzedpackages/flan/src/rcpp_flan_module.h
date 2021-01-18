#ifndef RCPP_FLAN_MODULE_H
#define RCPP_FLAN_MODULE_H

#include "FLAN_Sim.h"
#include "FLAN_MutationModel.h"
#include "FLAN_Clone.h"


RCPP_EXPOSED_AS(FLAN_Sim);
RCPP_EXPOSED_WRAP(FLAN_Sim);

RCPP_EXPOSED_AS(FLAN_SimClone);
RCPP_EXPOSED_WRAP(FLAN_SimClone);

RCPP_EXPOSED_AS(FLAN_SimInhomogeneous);
RCPP_EXPOSED_WRAP(FLAN_SimInhomogeneous);

RCPP_EXPOSED_AS(FLAN_ExponentialClone);
RCPP_EXPOSED_WRAP(FLAN_ExponentialClone);

RCPP_EXPOSED_AS(FLAN_DiracClone);
RCPP_EXPOSED_WRAP(FLAN_DiracClone);

RCPP_EXPOSED_AS(FLAN_InhomogeneousClone);
RCPP_EXPOSED_WRAP(FLAN_InhomogeneousClone);

RCPP_EXPOSED_AS(FLAN_MutationModel);
RCPP_EXPOSED_WRAP(FLAN_MutationModel);


#endif
