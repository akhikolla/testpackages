/*
 bcp: an R package for performing a Bayesian analysis
 of change point problems on a general graph.

 Copyright (C) 2012 Xiaofei Wang and Jay Emerson

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.r-project.org/Licenses/

 -------------------
 FILE: Cbcpgraph.cpp  */

/*  LIBRARIES  */
#include "bcp.h"

/* Define some structures/objects*/




void recomputeVals(Graph &graph, Partition &components, GraphParams &params)
{
  DoubleVec W(components.size(), 0.0);
  DoubleVec B(components.size(), 0.0);
  DoubleVec mean(params.kk, 0.0);
  DoubleMatrix means(components.size(), mean);
  int currblock, i, j;

  for (i = 0; i < graph.nodes.size(); i++) {
    currblock = graph.nodes[i].component;
    for (j = 0; j < params.kk; j++) {
      means[currblock][j] += graph.nodes[i].value[j];
      //   Rprintf("memb:%d, val:%0.2f means:%0.2f\n", currblock, graph.nodes[i].value, means[currblock]);
      W[currblock] += pow(graph.nodes[i].value[j], 2);
    }
  }

  for (i = 0; i < components.size(); i++) {
    B[i] = 0;

    for (j = 0; j < params.kk; j++) {
      means[i][j] /= components[i].size;
      B[i] += components[i].size * pow(means[i][j], 2);
    }

    W[i] -= B[i];
    Rprintf("Recomputed: i:%d, W: %0.2f, B: %0.2f, size: %d, %0.2f\n", i, W[i], B[i],
            components[i].size);
  }
}



void fullPixelPass(Graph &graph, Partition &components, GraphParams &params, MCMC &mcmc)
{
  int i, j, maxComp;
  int currblock, newblock;
  double maxll;

  for (i = 0; i < params.nn; i++) {
    // Rprintf("i:%d\n", i);
    currblock = graph.nodes[i].component;
    maxComp = components.size() + (components[currblock].size != graph.nodes[i].size);
    vector<Component> possibleBlocks = components; // this vector is not a partition!
    vector<MCMCStepGraph> possibleSteps(maxComp, mcmc.step);
    maxll = mcmc.step.lik;

    possibleBlocks[currblock].removeNode(graph.nodes[i]); // remove the node from current block

    // loop through all possible components
    for (j = 0; j < maxComp; j++) {
      if (j == components.size()) {
        Component newestBlock(graph.nodes[i]);
        possibleBlocks.push_back(newestBlock);
      } else if (j != currblock) {
        possibleBlocks[j].addNode(graph.nodes[i]);
      }

      possibleSteps[j].updateLogLik(params, graph, components, possibleBlocks[j],
                                    possibleBlocks[currblock], graph.nodes[i], j);

      if (possibleSteps[j].lik > maxll) {
        maxll = possibleSteps[j].lik;
      }
    }

    newblock = sampleLogLik(possibleSteps, maxll);
    updateComponents(params, mcmc, components, graph, possibleSteps, 
                     possibleBlocks, currblock,
                     newblock, i, -1);
  }
}

void activePixelPass(Graph &graph, Partition &components, GraphParams &params, 
                     MCMC &mcmc)
{
  int i, j, k, index;
  int currblock, newblock, passtype;
  double maxll, u;

  if (params.p1 == 1) {
    passtype = 1;  // correct APP
  } else if (params.p1 == 0) {
    passtype = 2;  // modified APP
  } else {
    u = Rf_runif(0.0, 1.0);

    if (u < params.p1) {
      passtype = 1;
    } else {
      passtype = 2;
    }
  }

  for (i = 0; i < params.nn; i++) {
    if (graph.nodes[i].active == 0) {
      continue;
    }

    //    Rprintf("node: %d  active: %d\n", i+1, graph.nodes[i].active );
    currblock = graph.nodes[i].component;
    vector<Component> possibleBlocks = components; // this vector is not a partition!
    vector<MCMCStepGraph> possibleSteps;
    IntegerVector neighbors = graph.nodes[i].neighbors;
    maxll = mcmc.step.lik;
    possibleBlocks[currblock].removeNode(graph.nodes[i]); // remove the node from current block

    IntVec allBlocks; // keep track of the components we are trying out (in the order we tried them)
    index = 0; // the unique count of the block we are on

    if (graph.nodes[i].active == 2) {
      mcmc.type2pix[mcmc.k - 101]++;
    }

    //if (graph.nodes[i].active == 1 || (graph.nodes[i].active == 2 && params.boundaryType == 1)) {
    if (graph.nodes[i].active == 1 || passtype == 2) {
      // if (graph.nodes[i].active == 1) {
      // agree with some neighbors, not all
      IntVec blockTried(possibleBlocks.size(), 0);

      // loop through all possible components
      for (k = 0; k < neighbors.size(); k++) {
        j = graph.nodes[neighbors[k]].component;

        if (blockTried[j] == 1) {
          continue;
        }

        if (j != currblock) {
          possibleBlocks[j].addNode(graph.nodes[i]);
        }

        allBlocks.push_back(j);
        blockTried[j] = 1;

        possibleSteps.push_back(mcmc.step);
        possibleSteps[index].updateLogLik(params, graph, components, possibleBlocks[j],
                                          possibleBlocks[currblock], graph.nodes[i], j);

        if (possibleSteps[index].lik > maxll) {
          maxll = possibleSteps[index].lik;
        }

        index++;
      }
    } else { // active = 2; disagrees with all neighbors
      IntVec blocksToTry(possibleBlocks.size(), 1); // indicators of whether or not we try these blocks

      for (k = 0; k < neighbors.size(); k++) {
        blocksToTry[graph.nodes[neighbors[k]].component] = 0;
        //Rprintf("nb: %d (%d)  ", neighbors[k]+1,  graph.nodes[neighbors[k]].component);
      }

      // loop through all possible components
      for (j = 0; j <= components.size(); j++) {
        if (j == components.size()) { // consider making own block
          if (components[currblock].size != graph.nodes[i].size) {
            Component newestBlock(graph.nodes[i]);
            possibleBlocks.push_back(newestBlock);
            allBlocks.push_back(j);
          } else {
            continue;
          }
        } else if (blocksToTry[j] == 0) {
          continue;
        } else {
          if (j != currblock) {
            possibleBlocks[j].addNode(graph.nodes[i]);
          }

          allBlocks.push_back(j);
        }

        possibleSteps.push_back(mcmc.step);
        possibleSteps[index].updateLogLik(params, graph, components, possibleBlocks[j],
                                          possibleBlocks[currblock], graph.nodes[i], j);

        if (possibleSteps[index].lik > maxll) {
          maxll = possibleSteps[index].lik;
        }

        index++;
      }
    }

    index = sampleLogLik(possibleSteps, maxll);

    newblock = allBlocks[index];
    //    Rprintf("currblock:%d, newblock:%d\n", currblock, newblock);
    updateComponents(params, mcmc, components, graph, possibleSteps, possibleBlocks, currblock,
                     newblock, i, index);
    //    if (graph.nodes[i].active == 2) Rprintf("i: %d  currblock: %d  newblock: %d   ll:%0.4f\n", i, currblock, newblock, mcmc.step.lik);

  }

}

// this computes the change in boundary length given a new block for pixel (i,k);
// returns 2x the length

// [[Rcpp::export]]
SEXP rcpp_ppm(SEXP pdata, SEXP pid, SEXP padj, SEXP pmcmcreturn, 
              SEXP pburnin, SEXP pmcmc, SEXP pa, SEXP pc,
              SEXP pmembs, SEXP pboundaryType, SEXP p1, SEXP pfreqAPP)
{

  NumericMatrix data(pdata);

  List adj(padj);
  NumericVector membInit(pmembs);

  int mcmcreturn = INTEGER_DATA(pmcmcreturn)[0];

  // some helper variables
  int i, m, j;
  double wstar;
  int freqAPP = INTEGER_DATA(pfreqAPP)[0];
  // initialize my graph
  // Rprintf("Initializing objects\n");
  Graph graph(pdata, pid, membInit, adj, false);
  GraphParams params(pc, pa, graph.nodes.size(), 
                data.nrow(), pboundaryType, pburnin,
                pmcmc, p1, data.ncol());
  int MM = params.burnin + params.mcmc;
  
  NumericMatrix pmean(params.nn, params.kk);
  NumericMatrix ss(params.nn, params.kk);
  NumericMatrix pvar(params.nn, params.kk);
  NumericVector pboundary(params.nn);
  // graph.print();
  // params.print();
  Partition components;
    
  // Rprintf("Initializing graph\n");

  for (i = 0; i < params.nn; i++) {
    graph.nodes[i].calcActiveAndBound(graph.nodes);

    if ((int) membInit[i] >= components.size()) {
      Component newComp(graph.nodes[i]);
      components.push_back(newComp);
    } else {
      components[(int) membInit[i]].addNode(graph.nodes[i]);
    }
  }

  MCMC mcmc(components, graph, params);
  // mcmc.step.print();
  // printPartition(components);
  /* --- MAKE THINGS TO BE RETURNED TO R --- */

  // Rprintf("Initializing matrix\n");
  int MM2, nn2;
  if (mcmcreturn == 1) {
    MM2 = MM;
    nn2 = params.nn;
  } else {
    MM2 = 1;
    nn2 = 1;
  }
  NumericMatrix membsa(nn2, MM2);
  NumericMatrix results(MM2*nn2,params.kk);
  double tmpMean;
  
  /* --- END THINGS TO BE RETURNED TO R --- */
  GetRNGstate(); // Consider Dirk's comment on this.
  // mcmc.step.print();
  // Rprintf("Before 100 pixel passes loop\n");
  /* ----------------START THE BIG LOOP--------------------------- */
  // first 100 iterations, do complete pixel passes
  for (m = 0; m < 100; m++) {
    fullPixelPass(graph, components, params, mcmc);
    mcmc.addStep(params);
    // mcmc.step.print();
  }

  // Rprintf("After 100 full pixel passes\n");
  for (m = 0; m < MM; m++) {
    fullPixelPass(graph, components, params, mcmc);
    for (i = 0; i < freqAPP; i++) {
      activePixelPass(graph, components, params, mcmc);
    }

    mcmc.addStep(params);
    // Rprintf("m:%d\n", m);
    // mcmc.step.print();
    // printPartition(components);
    wstar = mcmc.wstarvals[mcmc.k-1];  
    // Rprintf("1037 k:%d wstar:%0.8f\n", mcmc.k-1, wstar);
    if (m >= params.burnin || mcmcreturn == 1) {
      // for calculating posterior estimate of error variance
      // if (m >= burnin) {
      //   pvar += (mcmc.step.W + mcmc.step.B * wstar) / (params.nn2 * params.mm - 3);
      // }
      for (i = 0; i < params.nn; i++) {    
        if (m >= params.burnin)
          pboundary[i] += (graph.nodes[i].active > 0);
        if (mcmcreturn == 1)            
          membsa(i, m) = graph.nodes[i].component;
        for (j = 0; j < params.kk; j++) {
          tmpMean = (1 - wstar) * components[graph.nodes[i].component].mean[j]
                               + wstar * graph.mean;
                            
          if (m >= params.burnin) {
            pmean(i, j) += tmpMean;
            ss(i, j) += tmpMean * tmpMean;
          }
            
          if (mcmcreturn == 1) {
            results(m*params.nn+i, j) = tmpMean;
          }            
        }     
      }
    }
    

  }

  // mcmc.step.print();

  // Rprintf("Begin post-processing\n");

  for (i = 0; i < params.nn; i++) {
    pboundary[i] /= params.mcmc;
    for (j = 0; j < params.kk; j++) {
      pmean(i, j) /= params.mcmc;
      pvar(i, j) = (ss(i, j)/ params.mcmc - pmean(i,j)*pmean(i,j))*
        (params.mcmc/(params.mcmc-1));
    }
  }

  // pvar /= itsmcmc;

  PutRNGstate();
  List z;
  z["posterior.mean"] = pmean;
  z["posterior.var"] = pvar;
  // z["ll"] = wrap(mcmc.ll);
  z["posterior.prob"] = pboundary;
  z["mcmc.rhos"] = membsa;
  z["mcmc.means"] = results;
  z["blocks"] = wrap(mcmc.Mvals);
  z["len"] = wrap(mcmc.step.len);
  // z["type2pix"] = wrap(mcmc.type2pix);
  // z["simErr"] = wrap(mcmc.simErr);
  // Rprintf("End printing\n");

  //  closeLogFile();

  return z;
}
/* END MAIN  */

