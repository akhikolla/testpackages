#include "bcp.h"

using namespace Rcpp;

//simulate n samples from MVN(0, sigma=1e-5)
mat mvrnormArma(int n, Params& params) {
  mat Y = randn(n, params.kk);
  return Y * params.sigma_jitter;
}


double logKcalc(int bsize, int tau, Params& params) {
  double kratio = params.d/(bsize+params.d);
  int canFitReg = (bsize >= params.nreg)*1;
  double tmp = (kratio*canFitReg + 
                (1-canFitReg))*(tau==0) + 
                (1-kratio)*canFitReg*(tau==1);
  // Rprintf("tau:%d, bsize:%d K:%0.2f\n", tau, bsize, log(tmp));  
  return log(tmp);
}


double likelihood(double B, double W, int b, Params& params, double logC, 
                  double Q, double K,
                  int type) {
  // check if any w2 = 0 outside
  double lik;
  if (params.reg) { // type defaults to 1
    double Wtilde = W - Q;
    if (b == 1) {
      lik = logC + log(params.w[0])
      - (params.nn2 - 1) * log(Wtilde) / 2;
    } else {
      lik = logC - (b + 1) * log(B) / 2 - (params.nn2 - b - 2) * log(Wtilde) / 2
      + Rf_pbeta(B * params.w[0] / Wtilde / (1 + B * params.w[0] / Wtilde), (double) (b + 1) / 2,
                 (double) (params.nn2 - b - 2) / 2, 1, 1)
      + Rf_lbeta((double) (b + 1) / 2, (double) (params.nn2 - b - 2) / 2);
    }
    if (type == 1) {
      lik += K + params.priors[b-1];
    }
  } else {
    if (B == 0) {
      lik = params.priors[b - 1] + (params.kk * b + 1) * log(params.w[0]) / 2 -
        (params.nn2 * params.kk - 1) * log(W) / 2 - log((params.kk * b + 1)/2);
    } else if (b >= params.nn - 4 / params.kk) {
      lik = -DBL_MAX;
    } else {
      lik = params.priors[b - 1] - (params.kk * b + 1) * log(B) / 2 -
        ((params.nn2 - b) * params.kk - 2) * log(W) / 2
      + Rf_pbeta(B * params.w[0] / W / (1 + B * params.w[0] / W), 
                 (double) (params.kk * b + 1) / 2,
                 (double) ((params.nn2 - b) * params.kk - 2) / 2, 1, 1)
      + Rf_lbeta((double) (params.kk * b + 1) / 2,
                 (double) ((params.nn2 - b) * params.kk - 2) / 2);
    }
  }
  return lik;
}



// for graphs: binary search
int sampleLogLik(vector<MCMCStepGraph> possibleSteps, double maxll)
{
  int j;
  double myrand = Rf_runif(0.0, 1.0);
  DoubleVec llcum(possibleSteps.size());
  
  llcum[0] = exp(possibleSteps[0].lik - maxll);
  for (j = 1; j < possibleSteps.size(); j++) {
    llcum[j] = llcum[j - 1] + exp(possibleSteps[j].lik - maxll);
  }
  
  // begin binary search
  int startkk = 0;
  int endkk = llcum.size() - 1;
  int midkk;
  
  while (startkk != endkk) {
    midkk = floor((startkk + endkk) / 2);
    
    if (myrand <= llcum[midkk] / llcum[llcum.size() - 1]) {
      // search lower half
      endkk = midkk;
    } else {
      // search upper half
      startkk = midkk + 1;
    }
  }
  return (endkk);
}

void updateComponents(GraphParams &params, MCMC &mcmc, Partition &components, 
                      Graph &graph,
                      vector<MCMCStepGraph> &possibleSteps, 
                      vector<Component> &possibleBlocks, int currblock, int newblock,
                      int nodeId, int index = -1)
{
  int i;
  if ((params.reg && (newblock == currblock && 
        components[currblock].tau == possibleBlocks[index+1].tau)) || 
        (!params.reg && (newblock == currblock))) {
      return;
  } 
 if (index != -1) {
    mcmc.step = possibleSteps[index]; // this is if possiblesteps is a diff length
  } else { // multivariate only
    mcmc.step = possibleSteps[newblock]; 
  }

  if (params.reg) {
    if (newblock == currblock) {
      // tau changed
      components[currblock] = possibleBlocks[index+1];
      return;
    } else { // newly added 11.01 to check movedBlocks
      if (params.doneBurnin) mcmc.movedBlock[nodeId]++;
    }
    components[currblock] = possibleBlocks[0];
    if (newblock == components.size()) {
      components.push_back(possibleBlocks[index+1]);
    } else {
      components[newblock] = possibleBlocks[index+1];
    }
  } else {
    components[currblock] = possibleBlocks[currblock];
    
    if (newblock == components.size()) {
      components.push_back(possibleBlocks[newblock]);
    } else {
      components[newblock] = possibleBlocks[newblock];
    }
  }
  
  graph.updateNode(nodeId, newblock);
  
  // update the boundaryMatrix
  if (params.boundaryType == 1) {
    graph.updateBoundaryMatrix(nodeId, newblock, currblock);
  }
  
  if (components[currblock].size == 0) {
    if (currblock == components.size() - 1) {
      components.pop_back();
    } else {
      components[currblock] = components.back();
      components.pop_back();
      
      for (i = 0; i < graph.nodes.size(); i++) {
        if (graph.nodes[i].component == components.size()) {
          graph.nodes[i].component = currblock;
        }
        
        if (params.boundaryType == 1) {
          if (graph.boundarymat[components.size()][i] == 1) {
            graph.boundarymat[currblock][i] = 1;
            graph.boundarymat[components.size()][i] = 0;
          }
        }
      }
    }
  }
}

void printPartition(Partition &components)
{
  for (int i = 0; i < components.size(); i++) {
    Rprintf("i:%d ", i);
    components[i].print();
  }
}
