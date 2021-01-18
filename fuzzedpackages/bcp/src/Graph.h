#ifndef GRAPH_H
#define GRAPH_H


#include <vector>
#include <RcppArmadillo.h>  


/* DYNAMIC CREATION OF LOCAL VECTORS AND MATRICES */
typedef std::vector<double> DoubleVec;
typedef std::vector<int> IntVec;
typedef std::vector<IntVec> IntMatrix;
typedef std::vector<DoubleVec> DoubleMatrix;

using namespace std;
using namespace Rcpp;
using namespace arma;

// forward declarations (of functions actually in utils.cpp)
double logKcalc(int bsize, int tau, Params& params);

double likelihood(double B, double W, int b, Params& params, double logC=0, 
                  double Q=0, double K=-1.0, int type=1);
DoubleVec matrixCalcs(HelperVariables& helpers, 
                      Params& params, 
                      DoubleVec &w, 
                      int start, int end);
DoubleVec matrixCalcs(GraphParams& params, HelperVariables& helpers,
                      DoubleVec &w, 
                      uvec &nodes);
class Node
{
public:
  int id;
  DoubleVec value; // a vector, useful for multivariate
  int component;
  int active; // 0 = no, 1 = match some neighbor, 2 = disagree w all neighbors
  int boundlen;
  int size;
  IntegerVector neighbors;
  
  Node(DoubleVec &v, int c, int s, int idnum, List &adj)
  {
    id = idnum;
    value = v;
    component = c;
    active = 0;
    boundlen = 0;
    SEXP adjvec = adj[idnum];
    IntegerVector adj1(adjvec);
    neighbors = adj1;
    size = s;
  }

  // other methods
  void calcActiveAndBound(vector<Node> &nodes)
  {
    boundlen = 0;
    
    for (int i = 0; i < neighbors.size(); i++) {
      boundlen += (component != nodes[neighbors[i]].component);
    }
    
    if (boundlen == neighbors.size()) {
      active = 2;
    } else if (boundlen > 0) {
      active = 1;
    } else {
      active = 0;
    }
  }
  void printNeighbors(vector<Node> &nodes)
  {
    for (int i = 0; i < neighbors.size(); i++) {
      Rprintf("nb: %d, boundlen: %d, nb-comp:%d\n", neighbors[i], boundlen,
              nodes[neighbors[i]].component);
    }
  }
};


class Graph
{
public:
  vector<Node> nodes;
  double mean;
  IntMatrix boundarymat;
  uvec ids;
  double sumysq; // for W calculation
  
  
  Graph(SEXP &pdata, SEXP& pid, NumericVector &membs, List &adj, bool reg)
  {
    mean = 0.0;
    sumysq = 0.0;
    int M = 0;
    int i, j, nrow, nbComponent;
    
    // for multiple observations per location
    int nodesize = 0;
    int curr_id = 0;
    
    
    if (reg) {
      NumericVector data(pdata);
      ids = as<uvec>(pid);
      DoubleVec nodesum(1);
      nrow = data.size();
      
      for (i = 0; i < nrow; i++) {
        if (ids[i] > curr_id) { // on to next location
          Node node(nodesum, (int) membs[curr_id], nodesize, curr_id, adj);
          nodes.push_back(node);
          curr_id++;
          nodesum[0] = data[i];
          nodesize = 1;
          if (M < (int) membs[i] + 1) {
            M = (int) membs[i] + 1;
          }
        } else {
          nodesum[0] += data[i];
          nodesize++;
        }
        mean += data[i];    
        sumysq += pow(data[i], 2);
      }
      // the last node will need to be added
      Node node(nodesum, (int) membs[curr_id], nodesize, curr_id, adj);
      nodes.push_back(node);
      
      mean /= data.size();
      IntVec vboundary(nodes.size(), 0);
      boundarymat.assign(M, vboundary);
    } else {
      NumericMatrix data(pdata);
      DoubleVec values2;
      IntegerVector ids(pid);
      DoubleVec nodesum(data.ncol());
      
      for (i = 0; i < data.nrow(); i++) {
        if (ids[i] > curr_id) { // on to next location
          Node node(nodesum, (int) membs[curr_id], nodesize, curr_id, adj);
          nodes.push_back(node);
          curr_id++;
          nodesize = 1;
          if (M < (int) membs[i] + 1) {
            M = (int) membs[i] + 1;
          }
          for (j = 0; j < data.ncol(); j++) {
            mean += data(i, j);
            sumysq += pow(data(i, j), 2);
            nodesum[j] = data(i,j);
          }
        } else {
          nodesize++;
          for (j = 0; j < data.ncol(); j++) {
            mean += data(i, j);
            sumysq += pow(data(i, j), 2);
            nodesum[j] += data(i,j);
          }
        }

      }
      // the last node will need to be added
      Node node(nodesum, (int) membs[curr_id], nodesize, curr_id, adj);
      nodes.push_back(node);
      
      mean /= (data.nrow() * data.ncol());
      
      IntVec vboundary(data.nrow(), 0);
      boundarymat.assign(M, vboundary);
    }
    for (i = 0; i < nodes.size(); i++) {
      for (j = 0; j < nodes[i].neighbors.size(); j++) {
        nbComponent = nodes[nodes[i].neighbors[j]].component;
        
        if (nbComponent != nodes[i].component) {
          boundarymat[nbComponent][i] = 1;
        }
      }
    }
  }
  // other functions
  void print(bool all = false)
  {
    Rprintf("overall mean:%0.2f, overall ysq:%0.2f, num pixels: %d\n", 
            mean, sumysq, nodes.size());
    
    if (all) {
      for (int i = 0; i < nodes.size(); i++) {
        Rprintf("Node i:%d in block: %d, size:%d, sum(obs):%0.2f, boundlen: %d\n", i, nodes[i].component, 
                nodes[i].size, nodes[i].value[1], nodes[i].boundlen);
      }
      
      Rprintf("Boundary matrix\n");
      
      for (int i = 0; i < nodes.size(); i++) {
        for (int j = 0; j < 3; j++) {
          Rprintf("%d", boundarymat[j][i]);
        }
        
        Rprintf("\n");
      }
    }
  }
  void updateNode(int nodeId, int componentId)
  {
    nodes[nodeId].component = componentId;
    nodes[nodeId].calcActiveAndBound(nodes);
    
    for (int i = 0; i < nodes[nodeId].neighbors.size(); i++) {
      nodes[nodes[nodeId].neighbors[i]].calcActiveAndBound(nodes);
    }
  }
  void updateBoundaryMatrix(int nodeId, int newblock, int currblock)
  {
    boundarymat[newblock][nodeId] = 0;
    int boundaryOld = 0;
    int neighborNodeId, nbBoundaryOld, i, j;
    
    for (i = 0; i < nodes[nodeId].neighbors.size(); i++) {
      neighborNodeId = nodes[nodeId].neighbors[i];
      
      if (nodes[neighborNodeId].component == currblock) {
        boundaryOld = 1;
      }
      
      if (nodes[neighborNodeId].component != newblock) {
        boundarymat[newblock][neighborNodeId] = 1;
      }
      
      nbBoundaryOld = 0;
      
      for (j = 0; j < nodes[neighborNodeId].neighbors.size(); j++) {
        if (nodes[nodes[neighborNodeId].neighbors[j]].component == currblock
              && nodes[neighborNodeId].component != currblock) {
          nbBoundaryOld = 1;
          break;
        }
      }
      
      boundarymat[currblock][neighborNodeId] = nbBoundaryOld;
    }
    
    boundarymat[currblock][nodeId] = boundaryOld;
  }
  void recomputeBoundary(GraphParams &params, int M, int len) {
    // type =1 : node counting
    // type = 2: edge counting
    int blen = 0;
    int nbblock, i, j;
    
    if (params.boundaryType == 1) {
      IntVec vboundary(params.nn, 0);
      IntMatrix boundarymat2(M, vboundary);
      
      for (i = 0; i < nodes.size(); i++) {
        for (j = 0; j < nodes[i].neighbors.size(); j++) {
          nbblock = nodes[nodes[i].neighbors[j]].component;
          
          if (nodes[i].component != nbblock && boundarymat2[nbblock][i] == 0) {
            boundarymat2[nbblock][i] = 1;
            blen += 1;
          }
        }
      }
      
      for (i = 0; i < nodes.size(); i++) {
        for (j = 0; j < M; j++) {
          if (boundarymat2[j][i] != boundarymat[j][i]) {
            Rprintf("ERROR:\n");
          }
        }
      }
      
      if (blen != len) {
        Rprintf("ERROR len\n");
      }
      
    } else if (params.boundaryType == 2) {
      for (i = 0; i < nodes.size(); i++) {
        for (j = 0; j < nodes[i].neighbors.size(); j++) {
          blen += nodes[nodes[i].neighbors[j]].component != nodes[i].component;
        }
      }
    }
    // Rprintf("boundlen:%d\n", blen);
  }
  void checkBound(int M){// for error checking
    int totBound = 0;
    int totBound2 = 0;
    int i, j, nbBlock;
    for (i = 0; i < nodes.size(); i++) {
      IntVec blen(M, 0);
      
      for (j = 0; j < nodes[i].neighbors.size(); j++) {
        nbBlock = nodes[nodes[i].neighbors[j]].component;
        
        if (blen[nbBlock] == 0 && nbBlock != nodes[i].component) {
          blen[nbBlock] = 1;
          totBound++;
        }
      }
    }
    
    for (i = 0; i < nodes.size(); i++) {
      for (j = 0; j < M; j++) {
        totBound2 += boundarymat[j][i];
      }
    }
    
    Rprintf("totBound: %d | totBound2: %d\n", totBound, totBound2);
  }
};
class Component
{
  
public:
  int size;
  double Z;
  DoubleVec mean;
  double Q;
  double logC;
  double K;
  int tau;
  uvec nodeIds; // 1 indicates node is in block (length = N)
  uvec obsIds; // 1 indicates observation is in block (length = nn2)

  
  Component(GraphParams &params) { // multivariate constructor 
    size = 0;
    Z = 0.0;
    mean.assign(params.kk, 0);
  }
  Component(Node &node) { // multivariate constructor 
    size = node.size;
    Z = 0;
    
    for (int i = 0; i < node.value.size(); i++) {
      mean.push_back(node.value[i]/node.size);
      Z += pow(mean[i], 2);
    }
    Z *= size;
  }
  Component(GraphParams& params, Node &node, Graph &graph) { // regression
    size = node.size;
    mean = DoubleVec(1);
    mean[0] = node.value[0]/size;
    Z = pow(mean[0], 2);
    obsIds = zeros<uvec>(params.nn2);
    uvec these = find(graph.ids==node.id);
    for (int i = 0; i < these.n_rows; i++) {
      obsIds[these[i]] = 1;
    }
    nodeIds = zeros<uvec>(params.nn);
    nodeIds[node.id] = 1;
    tau = 0;
    Q = 0;
    K = logKcalc(size, tau, params);
    logC = 0;
  }
  
  // functions
  void initMemb(Node &node, Graph &graph)
  { // regression
    size += node.size;
    mean[0] = ((size - node.size) * mean[0] + node.value[0]) / size;
    Z = size * pow(mean[0], 2);
    uvec these = find(graph.ids==node.id);
    for (int i = 0; i < these.n_rows; i++) {
      obsIds[these[i]] = 1;
    }
    nodeIds[node.id] = 1;
  }
  void addNode(Node &node) { // multivariate
    size += node.size;
    Z = 0;
    for (int i = 0; i < node.value.size(); i++) {
      mean[i] = ((size - node.size) * mean[i] + node.value[i]) / size;
      Z += pow(mean[i], 2); 
    }
    Z *= size;
    // Rprintf("value:%0.2f, size:%0d, mean:%0.2f, Z:%0.2f\n", node.value[0],
    //   size, mean[0], Z);
  } 
  void addNode(GraphParams &params, HelperVariables &helpers, DoubleVec &w, 
            Node &node, Graph &graph, int ptau)
  { // regression
    size += node.size;
    mean[0] = ((size - node.size) * mean[0] + node.value[0]) / size;
    Z = size * pow(mean[0], 2);
    uvec these = find(graph.ids==node.id);
    for (int i = 0; i < these.n_rows; i++) {
      obsIds[these[i]] = 1;
    }
    nodeIds[node.id] = 1;
    changeTau(params, helpers, w, ptau);
  }
  void removeNode(Node &node) { // multivariate
    int i;
    // recompute the values for the current component (assuming we move this node to another component)
    Z = 0;
    if (size == node.size) {
      for (i = 0; i < node.value.size(); i++) {
        mean[i] = 0.0;
      }
      size = 0;
    } else {
      for (i = 0; i < node.value.size(); i++) {
        mean[i] = (mean[i] * size - node.value[i]) / (size - node.size);
        Z += pow(mean[i], 2);
      }
      size -= node.size;
      Z *= size;
    }
  }
  void removeNode(GraphParams &params, HelperVariables & helpers,
                  DoubleVec &w, Node &node, Graph &graph)
  { // regression
    // recompute the values for the current component (assuming we move this node to another component)
    if (size == node.size) {
      mean[0] = 0.0;
      Z = 0.0;
      size = 0;
    } else {
      mean[0] = (mean[0] * size - node.value[0]) / (size - node.size);
      size -= node.size;
      Z = size * pow(mean[0], 2);
    }
    uvec these = find(graph.ids==node.id);
    for (int i = 0; i < these.n_rows; i++) {
      obsIds[these[i]] = 0;
    }
    nodeIds[node.id] = 0;
    if (size < params.nreg) changeTau(params, helpers, w, 0);
    else changeTau(params, helpers, w, tau);
  }
  void changeTau(GraphParams &params, HelperVariables &helpers, 
                 DoubleVec &w, int ptau)
  {
    tau = ptau;
    K = logKcalc(size, tau, params);
    if (tau == 1) {
      DoubleVec out = matrixCalcs(params, helpers, w, obsIds);
      Q = out[0];
      logC = out[1]; 
    } else {
      Q = 0;
      logC = 0;
    }
  }
  void print()
  {
    Rprintf("Z: %0.2f, size:%d, mean: %0.2f Q:%0.2f logC:%0.2f K:%0.2f tau:%0d\n", 
            Z, size, mean[0], Q, logC, K, tau);
  }
};  

typedef vector<Component> Partition;



#endif
