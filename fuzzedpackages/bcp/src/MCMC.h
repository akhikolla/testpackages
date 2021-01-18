#ifndef MCMCSTEP_H
#define MCMCSTEP_H

#include <vector>
#include <RcppArmadillo.h>

#include "Graph.h"


/* DYNAMIC CREATION OF LOCAL VECTORS AND MATRICES */
typedef std::vector<double> DoubleVec;
typedef std::vector<int> IntVec;
typedef std::vector<IntVec> IntMatrix;
typedef std::vector<DoubleVec> DoubleMatrix;
static std::vector<double> DEFAULT_DOUBLEVEC;

using namespace std;
using namespace Rcpp;
using namespace arma;

// forward declarations (of functions actually in utils.cpp)
double logKcalc(int bsize, int tau, Params& params);

double likelihood(double B, double W, int b, Params& params, double logC,
                  double Q, double K, int type);


class MCMCStep {
public:
  double W;
  double B;
  int b;
  double lik;
  // constructors
  MCMCStep() {
    W = 0.0;
    B = 0.0;
    b = 1;
    lik = 0.0;
  }
  MCMCStep(const MCMCStep& step) {
    W = step.W;
    B = step.B;
    b = step.b;
    lik = step.lik;
  }

};

class MCMCStepGraph: public MCMCStep {
public:
  double K;
  double logC;
  double Q;

  DoubleVec w;
  int len;

  MCMCStepGraph() : MCMCStep() {
    len = 0;
    logC = 0.0;
    Q = 0.0;
    K = 0.0;
  }
  MCMCStepGraph(Partition &components, Graph &graph,
                           GraphParams &params,
                           DoubleVec &w0 = DEFAULT_DOUBLEVEC) : MCMCStep()
  {
    len = 0;
    Q = 0;
    logC = 0;
    K = 0;

    int i, j;
    W = graph.sumysq;
    B = -params.nn2 * pow(graph.mean, 2);
    if (!params.reg) B *= params.kk;
    b = components.size();

    for (i = 0; i < components.size(); i++) {
      W -= components[i].Z;
      B += components[i].Z;
      if (params.reg) {
        Q += components[i].Q;
        logC += components[i].logC;
        K += components[i].K;
      }

    }
    if (params.reg) w = w0;
    for (i = 0; i < params.nn; i++) {
      if (params.boundaryType == 1) {
        for (j = 0; j < b; j++) {
          len += graph.boundarymat[j][i];
        }
      }
      else if (params.boundaryType == 2)
        len += graph.nodes[i].boundlen;
    }
    // Rprintf("W:%0.2f, B: %0.2f\n", W, B);
    calcLogLik(params);
  }

  // other functions

  void calcLogLik(GraphParams &params)
  {
    if (abs(W) < 1.0e-12) {
      W = 1.0e-12;
    }

    if (params.reg) {
      double Wtilde = W - Q;

      if (b == 1) {
        lik = logC + K + log(params.w[0])
        - (params.nn2 - 1) * log(Wtilde) / 2;
      } else if (b >= params.nn - 5) {
        lik = -DBL_MAX;
      } else {
        double xmax = (B * params.w[0] / Wtilde) / (1 + (B * params.w[0] / Wtilde));
        lik = logC + K + len * log(params.p0)
          + Rf_pbeta(xmax, (double) (b + 1) / 2, (double) (params.nn2 - b - 2) / 2, 1, 1)
          + Rf_lbeta((double) (b + 1) / 2, (double) (params.nn2 - b - 2) / 2)
          - (b + 1) * log(B) / 2
          - (params.nn2 - b - 2) * log(Wtilde) / 2;
      }
    } else {
      double xmax = B * params.w[0] / W / (1 + B * params.w[0] / W);

      if (B == 0) {
        lik = len * log(params.p0) + (params.kk + 1) * log(params.w[0]) / 2
        - (params.kk * params.nn2 - 1) * log(W) / 2;
      } else if (b >= params.nn - 4 / params.kk) {
        lik = -DBL_MAX;
      } else {
        lik = len * log(params.p0)
        + Rf_pbeta(xmax, (double) (params.kk * b + 1) / 2,
                   (double) ((params.nn2 - b) * params.kk - 2) / 2, 1, 1)
        + Rf_lbeta((double) (params.kk * b + 1) / 2,
                   (double) ((params.nn2 - b) * params.kk - 2) / 2)
        - (params.kk * b + 1) * log(B) / 2
        - ((params.nn2 - b) * params.kk - 2) * log(W) / 2;
      }
    }

  }

  void updateLogLik(GraphParams &params, Graph &graph,
                    Partition &partition, Component &newcomp,
                    Component &oldcomp, Node &node, int newCompId)
  {
    if (params.reg) {
      int oldtau = partition[node.component].tau;
      if (newCompId == node.component && oldtau == newcomp.tau) {
        return;
      }
      if (newCompId != node.component) {
        int neighborsOldBlock, neighborOfOldComp, i, j;
        double Zdiff;
        b += 1 * (newCompId == b) - 1 * (partition[node.component].size == node.size);
        // update bound length
        if (params.boundaryType == 1) {
          if (newCompId >= graph.boundarymat.size()) {
            IntVec vboundary(params.nn, 0);
            graph.boundarymat.push_back(vboundary);
          }
          neighborsOldBlock = 0;
          for (i = 0; i < node.neighbors.size(); i++) {
            // subtract neighboring nodes that are no longer on the boundary of oldComp
            // this disqualifies any neighboring nodes that are in oldComp
            if (graph.nodes[node.neighbors[i]].component != node.component) {
              // these now used to be boundary of oldComp;
              // check if any of their other neighbors are in oldComp
              neighborOfOldComp = 0;
              Node neighborNode = graph.nodes[node.neighbors[i]];
              for (j = 0; j < neighborNode.neighbors.size(); j++) {
                if (neighborNode.neighbors[j] == node.id) {
                  continue;
                }
                if (graph.nodes[neighborNode.neighbors[j]].component == node.component) {
                  neighborOfOldComp = 1;
                  break;
                }
              }
              len -= (1 - neighborOfOldComp);
            } else {
              neighborsOldBlock = 1;
            }
            // add nodes that were not previously on the boundary of newComp but now are
            len -= graph.boundarymat[newCompId][node.neighbors[i]];
            len += (graph.nodes[node.neighbors[i]].component != newCompId);
          }
          // node is no longer boundary of newComp
          // add 1 if node is boundary of old block
          len -= graph.boundarymat[newCompId][node.id];
          len += neighborsOldBlock;

        } else if (params.boundaryType == 2) {
          len -= 2 * node.boundlen;
          for (i = 0; i < node.neighbors.size(); i++) {
            len += 2 * (graph.nodes[node.neighbors[i]].component != newCompId);
          }
        }

        // update globals
        Zdiff = partition[node.component].Z - newcomp.Z - oldcomp.Z;
        if (newCompId < partition.size()) {
          Zdiff += partition[newCompId].Z;
        }
        B -= Zdiff;
        W += Zdiff;
      }

      Q += newcomp.Q - partition[node.component].Q;
      K += newcomp.K - partition[node.component].K;
      logC += newcomp.logC - partition[node.component].logC;
      if (newCompId != node.component) {
        Q += oldcomp.Q;
        K += oldcomp.K;
        logC += oldcomp.logC;
        if (newCompId < partition.size()) {
          Q -= partition[newCompId].Q;
          K -= partition[newCompId].K;
          logC -= partition[newCompId].logC;
        }
      }
    } else { // multivariate
      if (newCompId == node.component) {
        return;
      }

      int neighborsOldBlock, neighborOfOldComp, i, j;
      double Zdiff;
      b += 1 * (newCompId == b) - 1 * (partition[node.component].size == node.size);

      // update bound length
      if (params.boundaryType == 1) {
        if (newCompId >= graph.boundarymat.size()) {
          IntVec vboundary(params.nn, 0);
          graph.boundarymat.push_back(vboundary);
        }

        neighborsOldBlock = 0;

        for (i = 0; i < node.neighbors.size(); i++) {
          // subtract neighboring nodes that are no longer on the boundary of oldComp
          // this disqualifies any neighboring nodes that are in oldComp
          if (graph.nodes[node.neighbors[i]].component != node.component) {
            // these now used to be boundary of oldComp;
            // check if any of their other neighbors are in oldComp
            neighborOfOldComp = 0;
            Node neighborNode = graph.nodes[node.neighbors[i]];

            for (j = 0; j < neighborNode.neighbors.size(); j++) {
              if (neighborNode.neighbors[j] == node.id) {
                continue;
              }

              if (graph.nodes[neighborNode.neighbors[j]].component == node.component) {
                neighborOfOldComp = 1;
                break;
              }
            }

            len -= (1 - neighborOfOldComp);
          } else {
            neighborsOldBlock = 1;
          }


          // add nodes that were not previously on the boundary of newComp but now are
          len -= graph.boundarymat[newCompId][node.neighbors[i]];
          len += (graph.nodes[node.neighbors[i]].component != newCompId);
        }

        // node is no longer boundary of newComp
        // add 1 if node is boundary of old block
        len -= graph.boundarymat[newCompId][node.id];
        len += neighborsOldBlock;

      } else if (params.boundaryType == 2) {
        len -= 2 * node.boundlen;

        for (i = 0; i < node.neighbors.size(); i++) {
          len += 2 * (graph.nodes[node.neighbors[i]].component != newCompId);
        }
      }

      Zdiff = partition[node.component].Z - newcomp.Z - oldcomp.Z;

      if (newCompId < partition.size()) {
        Zdiff += partition[newCompId].Z;
      }

      B -= Zdiff;
      W += Zdiff;

      if (params.kk == 1 && b == 1) B = 0;
    }
    calcLogLik(params);
  }

  void updateLogLikForMerge(GraphParams &params, Graph &graph,
                            Partition &partition, Component &newcomp,
                            int currblock, int newblock)
  {
    int i;
    double Zdiff;
    b--;
    // update bound length
    if (params.boundaryType == 1) {
      // subtract boundaries between the old and new comp (since merged now)
      for (i = 0; i < params.nn; i++) {
        if (newcomp.nodeIds[i] == 1) {
          len = len - graph.boundarymat[newblock][i] - graph.boundarymat[currblock][i];
        }
        // subtract 1 if a node is boundary to both blocks;
        if (graph.boundarymat[currblock][i] == 1 && graph.boundarymat[newblock][i] == 1) {
          len--;
        }

      }
    }

    // update globals
    Zdiff = partition[newblock].Z + partition[currblock].Z - newcomp.Z;
    B -= Zdiff;
    W += Zdiff;
    Q += newcomp.Q - partition[newblock].Q - partition[currblock].Q;
    K += newcomp.K - partition[newblock].K - partition[currblock].K;
    logC += newcomp.logC - partition[newblock].logC - partition[currblock].logC;
    calcLogLik(params);
  }
  void print() {
    Rprintf("lik:%0.2f, W:%0.2f, B:%0.2f, logC:%0.2f, K:%0.2f, Q:%0.2f, len =%d, b=%d\n",
            lik, W, B, logC, K, Q, len, b);
    for (int i = 0; i < w.size(); i++) Rprintf("w: %0.6f", w[i]);
    Rprintf("\n");
  }
};
class MCMCStepSeq: public MCMCStep {
public:
  double K;
  double logC;
  double Q;

  DoubleVec w;

  // blocks variables
  IntVec btau;
  IntVec rho;

  IntVec bend;
  IntVec bsize;
  DoubleVec bZ;
  DoubleVec blogC;
  DoubleVec bK;
  DoubleVec bQ;
  DoubleMatrix bmean;

  // constructors
  MCMCStepSeq(const MCMCStepSeq& step) : MCMCStep(step) {
    // W = step.W;
    // B = step.B;
    // b = step.b;
    logC = step.logC;
    w = step.w;
    K = step.K;
    // lik = step.lik;
    Q = step.Q;
  }

  MCMCStepSeq(HelperVariables &helpers, Params &params) : MCMCStep()
  {
    for (int i = 0; i < params.nn - 1; i++) {
      rho.push_back(0);
      if (params.reg) {
        if (i > 0 && i <= params.kk)
          w.push_back(params.w[i]/2);
        else if (i == 0)
          w.push_back(params.w[0]);
      }
    }
    rho.push_back(1);
    btau.push_back(0);
    bend.push_back(params.nn - 1);
    bsize.push_back(params.nn2);

    if (params.reg) {
      bZ.push_back(pow(helpers.cumy[params.nn - 1], 2) / params.nn2);
      bK.push_back(logKcalc(params.nn2, btau[0], params));
      bQ.push_back(0);
      blogC.push_back(0);
      W = helpers.cumysq[params.nn - 1] - bZ[0];
      logC = blogC[0];
      K = bK[0];
      Q = bQ[0];
      lik = likelihood(B, W, b, params, logC, Q, K);
    } else {
      double bZtmp = 0;
      DoubleVec bmean1(params.kk);
      for (int i = 0; i < params.kk; i++) {
        bmean1[i] = helpers.cumymat[i][params.nn - 1] / params.nn2;
        bZtmp += pow(bmean1[i], 2) * params.nn2;
      }
      bmean.assign(1, bmean1);
      bZ.push_back(bZtmp);
      B = bZtmp - params.nn2 * params.kk * pow(helpers.ybar, 2);
      W = helpers.ysqall - bZtmp;
      lik = likelihood(B, W, b, params);
    }
  }

  //other methods
  void print() {
    Rprintf("MCMCStep Info\n");
    Rprintf("B: %0.4f  W:%0.4f  b:%d   K: %0.2f  logC:%0.2f Q:%0.6f lik:%0.2f w:%0.8f\n", B, W,
            b, K, logC, Q, lik, w[1]);
    if (btau.size() == 0)
      return;
    for (int i = 0; i < btau.size(); i++) {
      Rprintf("i:%d   tau:%d  bend:%d  bsize:%d  bZ:%0.2f  bK:%0.2f bQ:%0.2f\n", i,
              btau[i], bend[i], bsize[i], bZ[i], bK[i], bQ[i]);
    }
  }
};

class MCMC
{ // only needed for graphs
public:
  MCMCStepGraph step;
  DoubleVec ll;
  IntVec Mvals;
  DoubleVec wstarvals;
  IntVec boundlens;
  DoubleVec simErr;
  IntVec type2pix;

  int k; // keeping track of which position we're at

  // posterior stuff (only needed for regression)
  vec pmeans;
  vec pvar;
  vec ss;
  DoubleVec pboundary;
  DoubleVec movedBlock; // to keep track # times a node moves

  // constructors
  MCMC(Partition &components, Graph &graph, GraphParams &params,
       DoubleVec &w0 = DEFAULT_DOUBLEVEC)
  {
    MCMCStepGraph step1(components, graph, params, w0);
    step = step1;

    if (params.reg) {
      pvar = zeros<vec>(params.nn);
      pmeans = zeros<vec>(params.nn);
      ss = zeros<vec>(params.nn);
      pboundary.assign(params.nn, 0);
      movedBlock.assign(params.nn, 0);
    }
    simErr.assign(params.nn, 0);

    int MM = params.burnin + params.mcmc + 101;
    ll.assign(MM, 0);
    Mvals.assign(MM, 0);
    wstarvals.assign(MM, 0);
    type2pix.assign(params.mcmc + params.burnin, 0);
    boundlens.assign(MM, 0);

    k = 0;
    addStep(params);
  }

  // other methods
  void addStep(GraphParams &params)
  {
    ll[k] = step.lik;
    Mvals[k] = step.b;
    boundlens[k] = step.len;

    double wstar = 0.0;
    // step.print();

    if (params.reg) {
      if (step.b > 1) {
        double Wtilde = step.W - step.Q;
        double xmax = step.B * params.w[0] / Wtilde /
        (1 + (step.B * params.w[0] / Wtilde));
        // Rprintf("W:%0.2f Q:%0.2f\n", step.W, step.Q);

        wstar = log(Wtilde) - log(step.B)
          + Rf_pbeta(xmax, (double) (step.b + 3) / 2, (double) (params.nn2 - step.b - 4) / 2, 1, 1)
          + Rf_lbeta((double) (step.b + 3) / 2, (double) (params.nn2 - step.b - 4) / 2)
          - Rf_pbeta(xmax, (double) (step.b + 1) / 2, (double) (params.nn2 - step.b - 2) / 2, 1, 1)
          - Rf_lbeta((double) (step.b + 1) / 2, (double) (params.nn2 - step.b - 2) / 2);
          wstar = exp(wstar);
      } else {
        wstar = params.w[0] / 2;
      }
    } else {
      double xmax = step.B * params.w[0] / step.W / (1 + (step.B * params.w[0] / step.W));
      if (step.B > 0) {
        //    Rprintf("xmax:%0.2f, log(W/B):%0.2f\n", xmax, log(step.W) - log(step.B));
        wstar = log(step.W) - log(step.B)
        + Rf_pbeta(xmax, (double) (step.b * params.kk + 3) / 2,
                   (double) ((params.nn2 - step.b) * params.kk - 4) / 2, 1, 1)
        + Rf_lbeta((double) (step.b * params.kk + 3) / 2,
                   (double) ((params.nn2 - step.b) * params.kk - 4) / 2)
        - Rf_pbeta(xmax, (double) (step.b * params.kk + 1) / 2,
                   (double) ((params.nn2 - step.b) * params.kk - 2) / 2, 1, 1)
        - Rf_lbeta((double) (step.b * params.kk + 1) / 2,
                   (double) ((params.nn2 - step.b) * params.kk - 2) / 2);
        wstar = exp(wstar);
      } else {
        wstar = params.w[0] * (step.b * params.kk + 1) / (step.b * params.kk + 3);
      }
    }

    wstarvals[k] = wstar;

    k++;
  }
  void postProcessing(GraphParams &params, int mcmc, mat & betaPosts)
  { // regression only
    for (int i = 0; i < params.nn; i++) {
      pmeans[i] /= mcmc;
      pboundary[i] /= mcmc;
      simErr[i] /= mcmc;
      movedBlock[i] /= mcmc*(params.freqAPP+1);
      pvar[i] = (ss[i] /mcmc - pmeans[i]*pmeans[i])*(mcmc/(mcmc-1));
    }
    betaPosts /= mcmc;
    betaPosts.cols(params.kk+1, betaPosts.n_cols-1) -=
      betaPosts.cols(0, params.kk)%betaPosts.cols(0, params.kk);
    //  Rprintf("mcmcits: %d\n", params.itsMCMC);
  }
};
#endif
