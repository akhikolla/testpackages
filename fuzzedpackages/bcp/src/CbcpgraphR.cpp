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
 FILE: CbcpgraphR.cpp  */
#include "bcp.h"

vec betaPostPredCalcs(GraphParams& params, HelperVariables& helpers, 
                      mat &XtildeTrue,
                      DoubleVec &w, 
                      uvec &these, uvec &these2) {

  int i;
  int n = these.size();
  mat Winv = zeros(params.kk, params.kk);
  mat Xtilde = XtildeTrue;
  mat XXW = Xtilde.t()*Xtilde;
  bool ok = FALSE;
  while(!ok) {
    for (i = 0; i < params.kk; i++) {
      if (XXW(i,i) < 1e-12) {
        Xtilde = XtildeTrue+ mvrnormArma(n, params);
        XXW = Xtilde.t()*Xtilde;
        break;
      }
      Winv(i,i) = XXW(i,i)*w[i+1]/(1-w[i+1]);
      if (i == params.kk-1) ok = TRUE;
    }
  }
  XXW = XXW + Winv;
  mat sumxy = Xtilde.t()*helpers.Y.elem(these); 
  vec betapost = XXW.i()*sumxy;
  return betapost;
}
DoubleVec matrixCalcs(GraphParams& params, HelperVariables& helpers,
                      DoubleVec &w, 
                      uvec &nodes) {
  int i;
  DoubleVec ret(2); // (Z, logC)
  uvec these = find(nodes==1);
  int n = these.size();
  mat Winv = zeros(params.kk, params.kk);
  mat Pmat = eye(n, n) - ones(n,n)/n; 
  mat Xtilde = Pmat*helpers.X.submat(these, helpers.pred_cols);
  mat XXW = Xtilde.t()*Xtilde;
  bool ok = FALSE;
  while(!ok) {
    for (i = 0; i < params.kk; i++) {
      if (XXW(i,i) < 1e-12) {
        Xtilde = Pmat*(helpers.X.submat(these, helpers.pred_cols)+ mvrnormArma(n, params));
        XXW = Xtilde.t()*Xtilde;
        break;
      }
      Winv(i,i) = XXW(i,i)*w[i+1]/(1-w[i+1]);
      if (i == params.kk-1) ok = TRUE;
    }
  }
  XXW = XXW + Winv;
  mat sumxy = Xtilde.t()*helpers.Y.elem(these);
  ret[0] = as_scalar(sumxy.t()*XXW.i()*sumxy);
  
  double detval, detsign;
  mat tmp = XXW*Winv.i();
  log_det(detval, detsign, tmp);

  ret[1] = -0.5*detval;
  return ret; 
}



void recomputeVals(Graph &graph, Partition &components)
{
  // DoubleVec W(components.size(), 0.0);
  DoubleVec B(components.size(), 0.0);
  DoubleVec means(components.size(), 0.0);
  int currblock, i;
  for (i = 0; i < graph.nodes.size(); i++) {
    currblock = graph.nodes[i].component;
    means[currblock] += graph.nodes[i].value[0];
    //   Rprintf("memb:%d, val:%0.2f means:%0.2f\n", currblock, graph.nodes[i].value, means[currblock]);
    // W[currblock] += graph.nodes[i].value * graph.nodes[i].value;
  }
  for (i = 0; i < components.size(); i++) {
    means[i] /= components[i].size;
    B[i] = components[i].size * pow(means[i], 2);
    // W[i] -= B[i];
    Rprintf("Recomputed: i:%d, B: %0.2f, size: %d, mean: %0.2f\n", i, B[i],
            components[i].size, means[i]);
  }
}

void updateComponentsForMerge(GraphParams &params, MCMC &mcmc, Partition &components, 
                      Graph &graph, MCMCStepGraph &possibleStep,
                      Component &possibleBlock, 
                      int currblock, int newblock)
{
  int i;
  if (newblock == currblock) {
    return;
  }
  mcmc.step = possibleStep;
  components[newblock] = possibleBlock;

  // update the boundaryMatrix
  if (params.boundaryType == 1) {
    for (i = 0; i < params.nn; i++) {
      if (components[newblock].nodeIds[i]==1) {
        graph.updateNode(i, newblock);
        graph.boundarymat[newblock][i] = 0;            
      } else if (graph.boundarymat[currblock][i] == 1) 
        graph.boundarymat[newblock][i] = 1;
      graph.boundarymat[currblock][i] = 0;   // currblock effectively nonexistent now
    }
  }
  if (currblock == components.size() - 1) {
    components.pop_back();
  } else {
    components[currblock] = components.back();
    components.pop_back();

    for (i = 0; i < params.nn; i++) {
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
  graph.recomputeBoundary(params, mcmc.step.b, mcmc.step.len);
}
void fullPixelPass(Graph &graph, Partition &components, GraphParams &params, 
                   MCMC &mcmc, HelperVariables &helpers,
                   bool silent)
{
  int i, j, maxComp, l, tau;
  int currblock, s;
  double maxll;
  for (i = 0; i< params.nn; i++) {  
    currblock = graph.nodes[i].component;
    maxComp = components.size() + (components[currblock].size != graph.nodes[i].size);
    IntVec blockNums;
    // this vector is not a partition, but rather stores the modified
    // candidate component information
    // [0] will always be the original component minus node i
    vector<Component> possibleBlocks; 
    vector<MCMCStepGraph> possibleSteps;
    maxll = mcmc.step.lik;

    possibleBlocks.push_back(components[currblock]);
    possibleBlocks[0].removeNode(params, helpers, mcmc.step.w, graph.nodes[i], 
                                 graph); // remove the node from current block
    l = 1;
    // loop through all possible components
    for (j = 0; j < maxComp; j++) {
      for (tau = 0; tau < 2; tau++) {
        if (j == components.size()){
          if (tau == 0) {
            Component newestBlock(params, graph.nodes[i], graph);
            possibleBlocks.push_back(newestBlock);
            possibleSteps.push_back(mcmc.step);
          } else break;
        } else if (j == currblock) {
          if (tau == 1 && components[currblock].size < params.nreg) {
            break;              
          }
          possibleBlocks.push_back(components[j]);
          possibleSteps.push_back(mcmc.step);
          if (tau != components[currblock].tau) {                
            possibleBlocks[l].changeTau(params, helpers, mcmc.step.w, tau);
          }      
        } else {
          if (tau == 0) {
            possibleSteps.push_back(mcmc.step);
            possibleBlocks.push_back(components[j]);
            possibleBlocks[l].addNode(params, helpers, mcmc.step.w, graph.nodes[i], graph, 0);
          } else if (possibleBlocks[l-1].size >= params.nreg) {
            possibleSteps.push_back(mcmc.step);
            possibleBlocks.push_back(possibleBlocks[l-1]);
            possibleBlocks[l].changeTau(params, helpers, mcmc.step.w, tau);
          } else break;
        }
        possibleSteps[l-1].updateLogLik(params, graph, components, 
                          possibleBlocks[l],
                          possibleBlocks[0], 
                          graph.nodes[i], j);
        blockNums.push_back(j);
        if (possibleSteps[l-1].lik > maxll) {
          maxll = possibleSteps[l-1].lik;
        }
        l++;
      
      }
    }

    s = sampleLogLik(possibleSteps, maxll);  
    updateComponents(params, mcmc, components, graph, possibleSteps, 
                         possibleBlocks, currblock, blockNums[s],
                         i, s);

  }
}
void blockMergePass(Graph &graph, Partition &components, GraphParams &params, 
                    HelperVariables &helpers,MCMC &mcmc, bool silent)
{
  int i, j, l;
  int s;
  double maxll;
  // Rprintf("len:%d\n", mcmc.step.len);
  for (i = 0; i< mcmc.step.b; i++) {

    vector<Component> possibleBlocks; 
    vector<MCMCStepGraph> possibleSteps;
    IntVec blockNums;
    maxll = mcmc.step.lik;
    l = 0;
    for (j = 0; j < mcmc.step.b; j++) {

      if (components[i].tau != components[j].tau) continue;
      // recalculate our block quantities
      possibleBlocks.push_back(components[i]);
      blockNums.push_back(j); // store the destination block
      possibleSteps.push_back(mcmc.step);
      if (i != j) {
        possibleBlocks[l].size = components[i].size + components[j].size;
        possibleBlocks[l].mean[0] = (components[i].mean[0]*components[i].size + components[j].mean[0]*components[j].size)/
                                  possibleBlocks[l].size;
        possibleBlocks[l].Z = possibleBlocks[l].size*pow(possibleBlocks[l].mean[0],2);
        possibleBlocks[l].tau = components[i].tau;
        possibleBlocks[l].nodeIds = components[i].nodeIds + components[j].nodeIds;
        possibleBlocks[l].obsIds = components[i].obsIds + components[j].obsIds;
        possibleBlocks[l].K = logKcalc(possibleBlocks[l].size, possibleBlocks[l].tau, params);
        if (components[j].tau==1) {
          DoubleVec out = matrixCalcs(params, helpers, mcmc.step.w, possibleBlocks[l].obsIds);
          possibleBlocks[l].Q = out[0];
          possibleBlocks[l].logC = out[1]; 
        } 
        possibleSteps[l].updateLogLikForMerge(params, graph, components, possibleBlocks[l], i, j);
        if (possibleSteps[l].lik > maxll) {
          maxll = possibleSteps[l].lik;
        }
      }    
      // Rprintf("j:%d  len:%d\n", j, possibleSteps[l].len);
      l++;

    }
    s = sampleLogLik(possibleSteps, maxll); 
    updateComponentsForMerge(params, mcmc, components, 
                      graph, possibleSteps[s],
                      possibleBlocks[s], 
                      i, blockNums[s]);
  }
}
void wPass(Partition &components, GraphParams &params, MCMC &mcmc, 
           HelperVariables &helpers)
{
  int i,j;
  double probs;
  for (i = 1; i< params.w.size(); i++) {
    vector<Component> possibleBlocks = components; 
    MCMCStepGraph candidateStep = mcmc.step;
    candidateStep.w = mcmc.step.w;
    candidateStep.w[i] += Rf_runif(-0.05*params.w[i], 0.05*params.w[i]);
    if (candidateStep.w[i] > params.w[i] || candidateStep.w[i] < 0) continue;
    candidateStep.Q = 0;
    candidateStep.logC = 0;

    for (j = 0; j < mcmc.step.b; j++) {
      possibleBlocks[j].changeTau(params, helpers, candidateStep.w, 
                                  possibleBlocks[j].tau); 
      candidateStep.Q += possibleBlocks[j].Q;
      candidateStep.logC += possibleBlocks[j].logC;
    }

    candidateStep.calcLogLik(params);
    probs = exp(candidateStep.lik-mcmc.step.lik);
    probs = probs/(1+probs);
    if (Rf_runif(0.0,1.0) < probs) {
      mcmc.step = candidateStep;
      components = possibleBlocks;
    } 
  }
}
void activePixelPass(Graph &graph, Partition &components, GraphParams &params, 
                     MCMC &mcmc, HelperVariables &helpers, bool silent)
{
  int i, j, k, l, s, tau;
  int currblock, passtype;
  double maxll, maxComp;

  if (params.p1 == 1) {
    passtype = 1;  // correct APP
  } else if (params.p1 == 0) {
    passtype = 2;  // modified APP
  } else {
    double u = Rf_runif(0.0, 1.0);
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
    // Rprintf("i:%d  ",i);

    currblock = graph.nodes[i].component;
    IntegerVector neighbors = graph.nodes[i].neighbors;
    maxll = mcmc.step.lik;

    // if (graph.nodes[i].active == 2) {
    //   mcmc.type2pix[mcmc.k - 101]++;
    // }
    vector<Component> possibleBlocks; 
    vector<MCMCStepGraph> possibleSteps; 
    IntVec blockNums;   
    possibleBlocks.push_back(components[currblock]);
    possibleBlocks[0].removeNode(params, helpers, mcmc.step.w, graph.nodes[i], graph); // remove the node from current block
    l = 1;

    if (graph.nodes[i].active == 1 || passtype == 2) {
      IntVec blockTried(mcmc.step.b, 0);
      // loop through all possible components
      for (k = 0; k < neighbors.size(); k++) {
        j = graph.nodes[neighbors[k]].component;
        if (blockTried[j] == 1) continue;
        blockTried[j] = 1;
        for (tau = 0; tau < 2; tau++) {
          if (j == currblock) {
            if (tau == 1 && components[currblock].size < params.nreg) {
              break;              
            }
            possibleBlocks.push_back(components[j]);
            possibleSteps.push_back(mcmc.step);
            if (tau != components[currblock].tau) {                
              possibleBlocks[l].changeTau(params, helpers, mcmc.step.w, tau);
            }      
          } else {
            if (tau == 0) {
              possibleSteps.push_back(mcmc.step);
              possibleBlocks.push_back(components[j]);
              possibleBlocks[l].addNode(params, helpers, mcmc.step.w, graph.nodes[i], graph, 0);
            } else if (possibleBlocks[l-1].size >= params.nreg) {
              possibleSteps.push_back(mcmc.step);
              possibleBlocks.push_back(possibleBlocks[l-1]);
              possibleBlocks[l].changeTau(params, helpers, mcmc.step.w, tau);
            } else break;
          }
          possibleSteps[l-1].updateLogLik(params, graph, components, 
                            possibleBlocks[l],
                            possibleBlocks[0], 
                            graph.nodes[i], j);
          blockNums.push_back(j);
          if (possibleSteps[l-1].lik > maxll) {
            maxll = possibleSteps[l-1].lik;
          }
          l++;
          // Rprintf("i:%d j:%d l:%d currblock:%d\n", i, j, l, currblock);
        } 
      }  
      
    } else { // active = 2; disagrees with all neighbors
      maxComp = (double) (mcmc.step.b+(components[currblock].size > graph.nodes[i].size));
      IntVec blocksToTry(maxComp, 1); // indicators of whether or not we try these blocks
      for (k = 0; k < neighbors.size(); k++) {
        blocksToTry[graph.nodes[neighbors[k]].component] = 0;
      }
      // loop through all possible components
      for (j = 0; j < maxComp; j++) {
        if (blocksToTry[j] == 0) {
          continue;
        }
        for (tau = 0; tau < 2; tau++) {
          if (j == components.size()){
            if (tau == 0) {
                Component newestBlock(params, graph.nodes[i], graph);
                possibleBlocks.push_back(newestBlock);
                possibleSteps.push_back(mcmc.step);
            } else break;
          } else if (j == currblock) {
            if (tau == 1 && components[currblock].size < params.nreg) {
              break;              
            }
            possibleBlocks.push_back(components[j]);
            possibleSteps.push_back(mcmc.step);
            if (tau != components[currblock].tau) {                
              possibleBlocks[l].changeTau(params, helpers, mcmc.step.w, tau);
            }      
          } else {
            
            if (tau == 0) {
              possibleSteps.push_back(mcmc.step);
              possibleBlocks.push_back(components[j]);
              possibleBlocks[l].addNode(params, helpers, mcmc.step.w, graph.nodes[i], graph, 0);
            } else if (possibleBlocks[l-1].size >= params.nreg) {
              possibleSteps.push_back(mcmc.step);
              possibleBlocks.push_back(possibleBlocks[l-1]);
              possibleBlocks[l].changeTau(params, helpers, mcmc.step.w, tau);
            } else break;
          }
          possibleSteps[l-1].updateLogLik(params, graph, components, 
                            possibleBlocks[l],
                            possibleBlocks[0], 
                            graph.nodes[i], j);
          blockNums.push_back(j);
          if (possibleSteps[l-1].lik > maxll) {
            maxll = possibleSteps[l-1].lik;
          }
          l++;
        // Rprintf("i:%d j:%d l:%d currblock:%d\n", i, j, l, currblock);
        }
      }
    }
    s = sampleLogLik(possibleSteps, maxll);
    updateComponents(params, mcmc, components, graph, possibleSteps, 
                     possibleBlocks, currblock, blockNums[s],
                     i, s);
  }

}

// [[Rcpp::export]]
SEXP rcpp_ppmR(SEXP py, SEXP px, SEXP pgrpinds, SEXP pid, SEXP padj, SEXP pmcmcreturn, 
              SEXP pburnin, SEXP pmcmc, SEXP pa, SEXP pc,
              SEXP pmembs, 
              SEXP pboundaryType, SEXP pba, 
              SEXP p1, SEXP pfreqAPP, SEXP pnreg)
{
  List adj(padj);
  NumericVector membInit(pmembs);

  int mcmcreturn = INTEGER_DATA(pmcmcreturn)[0];

  // some helper variables
  mat grpInds = as<mat>(pgrpinds);
  int i, m, bsize;
  double wstar;

  int freqAPP = INTEGER_DATA(pfreqAPP)[0];
  // initialize my graph
  // Rprintf("Initializing objects\n");
  Graph graph(py, pid, membInit, adj, true);
  HelperVariables helpers(py, px);
  GraphParams params(helpers.Y.n_rows, pc, pa, graph.nodes.size(), 
                pboundaryType, pburnin,
                pmcmc, p1, pfreqAPP, pba, pnreg);
  int MM = params.burnin + params.mcmc;
  
  // params.print();
  Partition components;
  // Rprintf("Initializing graph\n");
  for (i = 0; i < params.nn; i++) {
    graph.nodes[i].calcActiveAndBound(graph.nodes);
    if ((int) membInit[i] >= components.size()) {
      Component newComp(params, graph.nodes[i], graph);
      components.push_back(newComp);
    } else {
      components[(int) membInit[i]].initMemb(graph.nodes[i], graph);
    }
  }
  DoubleVec w0(params.kk+1);
  for (i = 0; i <= params.kk; i++) {
    if (i > 0 && i <= params.kk)
      w0[i] = params.w[i]/2;
    else if (i == 0) 
      w0[i] = params.w[0];
  }
  // Rprintf("ybar:%0.4f\n", graph.mean);
  for (i = 0; i < components.size(); i++) {
    components[i].changeTau(params, helpers, w0, 0);
    // components[i].print();
  }
  MCMC mcmc(components, graph, params, w0);
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
  NumericMatrix taus(nn2, MM2); 
  // store (fitted, intercept, slopes)
  mat results = zeros(MM2*nn2, 2+params.kk); // store conditional means

  mat betaposts = zeros(params.nn, (1+params.kk)*2);
  uvec interceptcols1(1);
  uvec interceptcols2(1);
  interceptcols1[0] = 0;
  interceptcols2[0] = params.kk+1;
  /* --- END THINGS TO BE RETURNED TO R --- */
  GetRNGstate(); // Consider Dirk's comment on this.

  // Rprintf("Before 100 pixel passes loop\n");
  // mcmc.step.print();

  /* ----------------START THE BIG LOOP--------------------------- */
  // first 100 iterations, do complete pixel passes
  bool silent = true;
  for (m = 0; m < 100; m++) {
    fullPixelPass(graph, components, params, mcmc, helpers, silent);
    blockMergePass(graph, components, params, helpers, mcmc, silent);
    wPass(components, params, mcmc, helpers);   
    mcmc.addStep(params);
  }
  // Rprintf("After 100 full pixel passes\n");

  // now do the rest:
  // 
  for (m = 0; m < MM; m++) {
    // Rprintf("m:%d\n", m);
    fullPixelPass(graph, components, params, mcmc, helpers, true);
    if (mcmc.step.len > 0) {
      for (i = 0; i < freqAPP; i++) {
        activePixelPass(graph, components, params, mcmc, helpers, true);
      }
      blockMergePass(graph, components, params, helpers, mcmc, silent);
    }
    wPass(components, params, mcmc, helpers);

    // passes are done, now store important info
    mcmc.addStep(params);
    // mcmc.step.print();
    if (m == params.burnin) params.doneBurnin = true; // silly but this is used in another func
    wstar = mcmc.wstarvals[mcmc.k-1];

    if (params.doneBurnin || mcmcreturn == 1) {
      // for posterior estimate of error variance
      // if (params.doneBurnin) 
      //   mcmc.pvar += (mcmc.step.W + mcmc.step.B * wstar) / (params.nn2 - 3);
      for (i = 0; i < params.nn; i++) {
        if (mcmcreturn == 1){
          membsa(i, m) = graph.nodes[i].component;        
          if (i < mcmc.step.b) {
            taus(i, m) = components[i].tau;
          }   
        }  
        if (params.doneBurnin) 
          mcmc.pboundary[i] += (graph.nodes[i].active > 0);
        if(i < mcmc.step.b) { 
          bsize = components[i].size;
          uvec these = find(components[i].obsIds == 1);
          uvec these2 = find(components[i].nodeIds == 1);
          uvec resultRows = these2 + params.nn*m;          
          mat grpIndMat = grpInds(these2, these);  
          vec intercepts = grpIndMat*(ones(bsize, 1)*(components[i].mean[0] * (1 - wstar) 
                           + graph.mean * wstar));
          vec fitted = intercepts;
          
          if (components[i].tau == 1) {
            mat Pmat = eye(bsize, bsize) - ones(bsize, bsize)/bsize; 
            mat Xtilde = Pmat*helpers.X.submat(these, helpers.pred_cols);
            vec betapost = betaPostPredCalcs(params, helpers, Xtilde, mcmc.step.w, these, 
                                        these2);
            mat betapostMat = repmat(betapost.t(), these2.n_rows, 1); 
            intercepts -= betapostMat*helpers.X.submat(these, helpers.pred_cols).t()*ones(bsize, 1)/bsize;
            fitted += grpIndMat*Xtilde*betapost;

            if (params.doneBurnin)  {              
              betaposts(these2, helpers.pred_cols) +=  betapostMat;
              betaposts(these2, helpers.pred_cols+params.kk+1) += betapostMat%betapostMat;
            }
            if (mcmcreturn == 1)
              results(resultRows, helpers.pred_cols+1) = betapostMat;
          } 
          if (mcmcreturn == 1) {
            results(resultRows, interceptcols1) = fitted;
            results(resultRows, interceptcols1+1) = intercepts;
          }
          if (params.doneBurnin) {
            mcmc.pmeans.elem(these2) += fitted;
            mcmc.ss.elem(these2) += fitted % fitted;
            betaposts(these2, interceptcols1) += intercepts;
            betaposts(these2, interceptcols2) += intercepts % intercepts;
          }
            
        } 
      }
      
    }

    

  }
  // mcmc.step.print();
  // printPartition(components);

  // Rprintf("Begin post-processing\n");

  mcmc.postProcessing(params, params.mcmc, betaposts);
  PutRNGstate();
  //
  List z;
  z["posterior.mean"] = wrap(mcmc.pmeans);
  z["posterior.var"] = wrap(mcmc.pvar);
 
  // z["ll"] = wrap(mcmc.ll);
  z["posterior.prob"] = wrap(mcmc.pboundary);
  z["mcmc.rhos"] = membsa;
  z["blocks"] = wrap(mcmc.Mvals);
  z["len"] = wrap(mcmc.step.len);
  // z["type2pix"] = wrap(mcmc.type2pix);
  // z["simErr"] = wrap(mcmc.simErr);
  // z["tau"] = wrap(taus);
  z["mcmc.means"] = wrap(results);
  z["betaposts"] = wrap(betaposts);
  z["movedBlock"] = wrap(mcmc.movedBlock);
  // Rprintf("End printing\n");

  //  closeLogFile();

  return z;
}
/* END MAIN  */

