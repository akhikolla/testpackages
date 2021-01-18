#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector pruneEdgesHelper(IntegerVector support,
                               IntegerMatrix edgeList,
                               int numLeaves,
                               int numVertices) {
  IntegerVector supportToRet(support.size());
  IntegerMatrix edgeListToRet(edgeList.nrow(), 2);
  for(int i = 0; i < supportToRet.size(); i++) {
    supportToRet(i) = 1 * support(i);
  }
  for(int i = 0; i < edgeList.nrow(); i++) {
    for(int j = 0; j < 2; j++) {
      edgeListToRet(i, j) = 1.0 * edgeList(i, j) - 1;
    }
  }
  int numEdgesInSupport = 0;
  for (int i = 0; i < supportToRet.size(); i++) {
    numEdgesInSupport += supportToRet(i);
  }
  int k = 0;
  IntegerMatrix reducedEdgeList(numEdgesInSupport, 2);
  for (int i = 0; i < supportToRet.size(); i++) {
    if (supportToRet(i) == 1.0) {
      reducedEdgeList(k, 0) = edgeListToRet(i, 0);
      reducedEdgeList(k, 1) = edgeListToRet(i, 1);
      k++;
    }
  }
  IntegerVector vertexDegrees = rep(0, numVertices);
  for(int i = 0; i < reducedEdgeList.nrow(); i++) {
    for(int j = 0; j < reducedEdgeList.ncol(); j++) {
      vertexDegrees(reducedEdgeList(i, j))++;
    }
  }

  IntegerMatrix edgeToEdgeNumberMatrix(numVertices, numVertices);
  for(int i = 0; i < edgeListToRet.nrow(); i++) {
    edgeToEdgeNumberMatrix(edgeListToRet(i, 0), edgeListToRet(i, 1)) = i;
    edgeToEdgeNumberMatrix(edgeListToRet(i, 1), edgeListToRet(i, 0)) = i;
  }
  IntegerMatrix reducedAdjList(numVertices, 3);
  for(int i = 0; i < numVertices; i++) {
    for(int j = 0; j < 3; j++) {
      reducedAdjList(i, j) = -1;
    }
  }
  for(int i = 0; i < reducedEdgeList.nrow(); i++) {
    int a = reducedEdgeList(i, 0);
    int b = reducedEdgeList(i, 1);

    if (reducedAdjList(a, 0) == -1) { reducedAdjList(a, 0) = b; }
    else if(reducedAdjList(a, 1) == -1) { reducedAdjList(a, 1) = b; }
    else { reducedAdjList(a, 2) = b; }

    if (reducedAdjList(b, 0) == -1) { reducedAdjList(b, 0) = a; }
    else if(reducedAdjList(b, 1) == -1) { reducedAdjList(b, 1) = a; }
    else { reducedAdjList(b, 2) = a; }
  }

  int numDegOneInternalNodes = 0;
  IntegerVector degreeOneInternalNodes = rep(0, 2*(numVertices - numLeaves));
  for(int i = numLeaves; i < numVertices; i++) {
    if(vertexDegrees(i) == 1.0) {
      degreeOneInternalNodes(numDegOneInternalNodes) = i;
      numDegOneInternalNodes++;
    }
  }

  while(numDegOneInternalNodes != 0) {
    int curNode = degreeOneInternalNodes(numDegOneInternalNodes - 1);
    numDegOneInternalNodes--;
    if(vertexDegrees(curNode) == 0) {
      continue;
    }
    int neighbor = std::max(std::max(reducedAdjList(curNode, 0), reducedAdjList(curNode, 1)), reducedAdjList(curNode, 2));

    vertexDegrees(curNode)--;
    vertexDegrees(neighbor)--;

    reducedAdjList(curNode, 0) = -1;
    reducedAdjList(curNode, 1) = -1;
    reducedAdjList(curNode, 2) = -1;

    if(reducedAdjList(neighbor, 0) == curNode) { reducedAdjList(neighbor, 0) = -1; }
    else if(reducedAdjList(neighbor, 1) == curNode) { reducedAdjList(neighbor, 1) = -1; }
    else{ reducedAdjList(neighbor, 2) = -1; }

    supportToRet(edgeToEdgeNumberMatrix(curNode, neighbor)) = 0;

    if (vertexDegrees(neighbor) == 1 && neighbor >= numLeaves) {
      numDegOneInternalNodes++;
      degreeOneInternalNodes(numDegOneInternalNodes - 1) = neighbor;
    }
  }
  return supportToRet;
}
