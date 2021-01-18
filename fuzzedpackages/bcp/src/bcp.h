#ifndef BCP_H
#define BCP_H

/*  LIBRARIES  */
#include <RcppArmadillo.h>  
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R_ext/Random.h>
#include <R.h>
#include <Rdefines.h>
#include <vector>

/* DYNAMIC CREATION OF LOCAL VECTORS AND MATRICES */
typedef std::vector<double> DoubleVec;
typedef std::vector<int> IntVec;
typedef std::vector<DoubleVec> DoubleMatrix;
typedef std::vector<IntVec> IntMatrix;


#include "HelperVariables.h"
#include "Params.h"
#include "MCMC.h"


// =======================================  //
// ===  Declaring Functions  ===  //
// =======================================  //

// from utils.cpp (shared functions)

mat mvrnormArma(int n, Params& params);

double logKcalc(int bsize, int tau, Params& params);

double likelihood(double B, double W, int b, Params& params, double logC, 
                  double Q, double K,
                  int type);

int sampleLogLik(vector<MCMCStepGraph> possibleSteps, double maxll);

// int sampleFromLikelihoods(DoubleVec &likvals, double maxlik);

void printPartition(Partition &components);

void updateComponents(GraphParams &params, MCMC &mcmc, Partition &components, 
                      Graph &graph,
                      vector<MCMCStepGraph> &possibleSteps, 
                      vector<Component> &possibleBlocks, int currblock, int newblock,
                      int nodeId, int index);

// I will comment out the functions that must not be shared...

// from Cbcp.cpp

DoubleVec matrixCalcs(HelperVariables& helpers, 
                      Params& params, 
                      DoubleVec &w, 
                      int start, int end);

// MCMCStepSeq pass(MCMCStepSeq &step, HelperVariables &helpers, 
//               Params &params, int printtmp);

// from CbcpM.cpp
int sampleFromLikelihoodsM(DoubleVec &likvals);
// MCMCStepSeq pass(MCMCStepSeq &step, HelperVariables &helpers, Params &params);



// from CbcpgraphR.cpp

vec betaPostPredCalcs(GraphParams& params, HelperVariables& helpers, 
                      mat &XtildeTrue,
                      DoubleVec &w, 
                      uvec &these, uvec &these2);
DoubleVec matrixCalcs(GraphParams& params, HelperVariables& helpers,
                      DoubleVec &w, 
                      uvec &nodes);                     
void recomputeVals(Graph &graph, Partition &components);

void updateComponentsForMerge(GraphParams &params, MCMC &mcmc, Partition &components, 
                      Graph &graph, MCMCStepGraph &possibleStep,
                      Component &possibleBlock, 
                      int currblock, int newblock);
// void fullPixelPass(Graph &graph, Partition &components, GraphParams &params, 
//                    MCMC &mcmc, HelperVariables &helpers,
//                    bool silent);

// void blockMergePass(Graph &graph, Partition &components, GraphParams &params, 
//                     HelperVariables &helpers,MCMC &mcmc, bool silent);     

// void wPass(Partition &components, GraphParams &params, MCMC &mcmc, HelperVariables &helpers);

// void activePixelPass(Graph &graph, Partition &components, GraphParams &params, 
//                      MCMC &mcmc, HelperVariables &helpers, bool silent);


// from Cbcpgraph.cpp

void recomputeVals(Graph &graph, Partition &components, GraphParams &params);
// void fullPixelPass(Graph &graph, Partition &components, GraphParams &params, MCMC &mcmc);
// void activePixelPass(Graph &graph, Partition &components, GraphParams &params, 
//                      MCMC &mcmc);

#endif
