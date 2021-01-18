//  EFT_DQP.cpp
//    Exact functional test implementaiton using dynamic programming
//
//  Created by Hien Nguyen on 7/24/2018.
//
//  Revision history:
//  2019-09-07 (Hien Nguyen): fix the error case when the branch and bound occurs at the root node.
//
//  2019-02-25 (Hien Nguyen): insulate code into the namespace DP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.cpp of 2.4.7 to EFT_DQP.cpp to distinguish from other
//     implementations of the exact functional test.
//

//#include "Node.h"
#include "EFT_DQP.h"
#include "EFT_QP.h"
#include "trimTable.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>
#include <algorithm>
#include <climits>

using namespace std;

//convertToInt: convert Rs to equivalent integer
int DQP::convertToInt(vector<int> Rs) {
  std::sort(Rs.begin(), Rs.end());

  int eInt = 0;
  for (size_t x = 0; x < Rs.size(); x++) {
      eInt *= 127;   // MS: TODO 127 shall be changed to the maximum column sum. 9/22/2019
    eInt += Rs[x];
  }
  return eInt;
}


// isMember: to find the node in a layer
// if the node exists in the layer, return its index
// otherwise, return -1

int DQP::isMember(int &eInt, vector<Node> &layer) {
  for (size_t x = 0; x < layer.size(); x++) {
    if (eInt == layer[x].getEquiInt()) return (int) x;
  }
  return -1;
}

int DQP::searchHashTable(vector<vector< pair<int, int>>> &hashTable, int element) {
  int i = element % hashTable.size();
  for (size_t j = 0; j < hashTable[i].size(); j++) {
    if (hashTable[i][j].first == element) return hashTable[i][j].second;
  }
  return -1;
}

// compute the weight between two nodes
double DQP::colChisq(vector<int> &Rs1, vector<int> &Rs2, int sum, const vector<int> & squares, double &COLMARGIN) {
  if (sum  > 0) {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs2.size(); x++) {
      colchisq += squares[Rs1[x] - Rs2[x]];
    }
    colchisq =  colchisq * COLMARGIN / sum;
    return (colchisq);
  }
  else return 0;
}

double DQP::colChisq(vector<int> Rs1, int &sum, const vector<int> & squares, double &COLMARGIN) {
  if (sum  > 0) {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs1.size(); x++) {
      colchisq += squares[Rs1[x]];
    }
    colchisq =  colchisq * COLMARGIN / sum;
    return (colchisq);
  }
  else return 0;

}
// compute the length from the current node to the end node
double DQP::length(vector<int> &Rs1, int &sum, int &layer, vector<int> &Cs, const vector<double> & factorials) {
  double length = factorials[sum];
  for (size_t x = 0; x < Rs1.size(); x++) {
    length /= factorials[Rs1[x]];
  }
  for (int x = 0; x < layer; x++) {
    length /= factorials[Cs[x]];
  }
  return (length);
}

// compute the length between two nodes
double DQP::length(vector<int> &Rs1, vector<int> &Rs2, const vector<double> & factorials) {
  double length = 1.0;
  for (size_t x = 0; x < Rs2.size(); x++) {
    length /= factorials[Rs1[x]- Rs2[x]];
  }
  return (length);
}

// compute the funchisq without the fixed marginals
double DQP::funchisqByCol(vector<vector<int>> &observedTable,
                          vector<int> &CSUM, const vector<int> & squares, double &COLMARGIN) {
  double funchisq = 0.0;
  double colchisq = 0.0;
  for (size_t j = 0; j < observedTable[0].size(); j++) {
    colchisq = 0;
    if (CSUM[j] > 0) {
    for (size_t i = 0; i < observedTable.size(); i++) {
      colchisq += squares[observedTable[i][j]];
    }
    colchisq = colchisq * COLMARGIN / CSUM[j];
    }

    funchisq += colchisq;
  }
  return (funchisq);
}


double DQP::lower_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums, double &COLMARGIN)
{
  double lower_bound = 0;

  vector<int> U(Rsum);

  size_t nrows = Rsum.size();

  vector<size_t> order(nrows);
  for (size_t q = 0; q < nrows; ++q) {
    order[q] = q;
  }

  // sort U in increasing order
  sort(order.begin(), order.end(),
       [&U](size_t i1, size_t i2) {return U[i1] < U[i2]; });

  for (int l = 0; l < layer; l++) {
    // find lower bound for row l
    int runsum = 0;

    for (size_t k = 0; k < nrows; ++k) {
      // accumulate the lower bound
      double xavg = (O_colsums[l] - runsum) / (double)(nrows - k);
      if (U[order[k]] < xavg) {
        if  (O_colsums[l] > 0) lower_bound += U[order[k]] * U[order[k]] * COLMARGIN / (double)O_colsums[l];
        runsum += U[order[k]];
      }
      else {
        if  (O_colsums[l] > 0) lower_bound += (nrows - k)  * xavg  * xavg  * COLMARGIN / (double)O_colsums[l];
        break;
      }
    }
  }

  return lower_bound;
}


double DQP::upper_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums, double &COLMARGIN)
{
  double upper_bound = 0;

  vector<int> U(Rsum);

  size_t nrows = Rsum.size();

  vector<size_t> order(nrows);
  for (size_t q = 0; q < nrows; ++q) {
    order[q] = q;
  }

  // sort U in decreasing order
  sort(order.begin(), order.end(),
       [&U](size_t i1, size_t i2) {return U[i1] > U[i2]; });

  for (int l = layer - 1; l >= 0; l--) {
    // find lower bound for row l
    if (O_colsums[l] > 0) {
      int runsum = 0;

      for (size_t k = 0; k < nrows; ++k) {
        // accumulate the lower bound
        int xmax = O_colsums[l] - runsum;
        if (U[order[k]] < xmax) {
          if  (O_colsums[l] > 0) upper_bound += U[order[k]] * U[order[k]] * COLMARGIN / O_colsums[l];
          runsum += U[order[k]];
        }
        else if (xmax != 0) {
          if  (O_colsums[l] > 0)  upper_bound += xmax * xmax * COLMARGIN / O_colsums[l];
          runsum += xmax;
        }
        else {
          break;
        }
      }
    }
  }

  return upper_bound;

}


// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void DQP::createNode(Node &node, vector<int> Rs, vector<int> &Cs, int layer, vector<int> &currRs, int &nrows, int sum1, int sum2,
                     vector<int> &S, const int &i, const vector<int> & squares, const vector<double> & factorials, vector<Node> &Layer,
                     double &COLMARGIN, vector<vector<pair<int, int>>> &hashTable) {
  if (i == nrows) {
    double len = DQP::length(Rs, currRs, factorials);
    int colchisq = DQP::colChisq(Rs, currRs, Cs[layer], squares, COLMARGIN);

    int eInt = DQP::convertToInt(currRs);
    int index = DQP::searchHashTable(hashTable, eInt);

    // if the child node does not exist yet, insert it to the next layer as a new node
    if (index < 0) {
      Layer.push_back(Node(currRs, eInt));
      node.addChildLink((int) Layer.size() - 1, len, colchisq);
      //update hashTable by insertion
      hashTable[eInt % hashTable.size()].push_back(make_pair(eInt,Layer.size() - 1));

      Layer[Layer.size() - 1].setMinPastChisq(node.getMinPastChisq() + colchisq);
      Layer[Layer.size() - 1].setMaxPastChisq(node.getMaxPastChisq() + colchisq);

      Layer[Layer.size() - 1].setLB(DQP::lower_bound(layer, currRs, Cs, COLMARGIN));
      Layer[Layer.size() - 1].setUB(DQP::upper_bound(layer, currRs, Cs, COLMARGIN));

      Layer[Layer.size() - 1].setLengthToEnd(length(currRs, S[layer-1], layer, Cs, factorials));

    }

    // if the child node already exists, add a new link from the current node to that child node
    else {
      node.addChildLink(index, len, colchisq);

      Layer[index].setMinPastChisq(std::min(Layer[index].getMinPastChisq(), node.getMinPastChisq() + colchisq));
      Layer[index].setMaxPastChisq(std::max(Layer[index].getMaxPastChisq(), node.getMaxPastChisq() + colchisq));
    }

  }
  else {
    int lowerbound, upperbound;

    sum1 += (i > 0 ? Rs[i - 1] : 0);
    sum2 += (i > 0 ? currRs[i - 1]: 0);

    lowerbound = std::max(0, Rs[i]- Cs[layer] + sum1 - sum2);
    upperbound = std::min(Rs[i], (layer-1 >= 0 ? S[layer-1]: 0) - sum2);

    for (int x = lowerbound; x <= upperbound; x++) {
      currRs[i] = x;
      DQP::createNode(node, Rs, Cs, layer, currRs, nrows, sum1, sum2, S, i + 1, squares, factorials, Layer, COLMARGIN, hashTable);
    }
  }
}

double DQP::EFTNetwork(vector<vector<int>> observedTable)
{
  observedTable = trimTable(observedTable);

  int nrows = (int) observedTable.size();
  int ncols = nrows > 0 ? (int) observedTable[0].size() : 0;

  if(nrows == 0 && ncols == 0) {
    return 1.0;
  }

  int N = 0;
  vector<int> GlobalRowSums(nrows, 0);
  vector<int> ColSums(ncols, 0);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      N += observedTable[i][j];
      GlobalRowSums[i] += observedTable[i][j];
      ColSums[j] += observedTable[i][j];
    }
  }

  vector<int> squares(N);
  for (int x = 0; x < N; x++) squares[x] = x*x;

  vector<double> factorials(N+1);
  factorials[0] = 1.0;
  for (int x = 1; x <= N; x++) factorials[x] = x*factorials[x - 1];
  //for (int x = 1; x <= N; x++)  factorials[x] = factorial<double>(x);

  double marginal = factorials[N];
  for (int x = 0; x < nrows; x++) {
    marginal /= factorials[GlobalRowSums[x]];
  }
  for (int x = 0; x < ncols; x++) {
    marginal /= factorials[ColSums[x]];
  }

  std::vector<int> S(ncols);
  S[0] = ColSums[0];
  for (int x = 1; x < ncols; x++) {
    S[x] = S[x - 1] + ColSums[x];
  }

  double COLMARGIN = 1;
  for (int i = 0; i < ncols; i++) {
    if (ColSums[i]  > 0 ) COLMARGIN *= ColSums[i];
  }

  // compute the adjusted funchisq
  double funchisq = DQP::funchisqByCol(observedTable, ColSums, squares,COLMARGIN);

  vector<vector<Node>> network(ncols + 1);

  //create the start node
  network[ncols].push_back(Node(GlobalRowSums, 0));
  network[ncols][0].addPastLen(1.0, 0);
  network[ncols][0].setMaxPastChisq(0);
  network[ncols][0].setMinPastChisq(0);

  network[ncols][0].setLB(DQP::lower_bound(ncols, GlobalRowSums, ColSums, COLMARGIN));
  network[ncols][0].setUB(DQP::upper_bound(ncols, GlobalRowSums, ColSums, COLMARGIN));
  network[ncols][0].setLengthToEnd(DQP::length(GlobalRowSums, S[ncols-1], ncols, ColSums, factorials));

  if (network[ncols][0].getLB() >= funchisq)  return 1;

  // hash table
  int hashTableSize = 199;
  vector<vector< pair<int, int>>> hashTable(hashTableSize);

  // generate the network
  vector<int> currRs(nrows);
  for (int layer = ncols; layer > 1; layer--) {
    for (size_t n = 0; n < network[layer].size(); n++) {
      if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() < funchisq) {
        if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() >= funchisq) {
          DQP::createNode(network[layer][n], network[layer][n].getRsum(), ColSums, layer-1, currRs, nrows, 0, 0, S, 0, squares, factorials, network[layer - 1], COLMARGIN, hashTable);
        }
      }
    }

    // for (int x = 0; x < hashTableSize; x++) hashTable[x].clear(); //clear hashtable
  }

  double colchisq;
  for (size_t x = 0; x < network[1].size(); x++) {
    colchisq = colChisq(network[1][x].getRsum(), ColSums[0], squares, COLMARGIN);
    network[1][x].setLB(colchisq);
    network[1][x].setUB(colchisq);
  }

  // compute upperbound and lowerbound for higher layer
  double minLB = 0;
  double maxUB = 0;
  double tempBound;

  for (int layer = 2; layer <= ncols; layer++) {
    for (size_t node = 0; node < network[layer].size(); node++) {

      //network[layer][node].setLengthToEnd(length(network[layer][node].getRsum(), S[layer - 1], layer, ColSums, factorials));

      minLB = INT_MAX; maxUB = 0;

      if (network[layer][node].getSize() > 0) {

        for (int child = 0; child < network[layer][node].getSize(); child++) {

          colchisq = network[layer][node].getColChisqToChildren(child);

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getLB();
          if (minLB > tempBound) minLB = tempBound;

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getUB();
          if (maxUB < tempBound) maxUB = tempBound;

        }

        network[layer][node].setLB(minLB);
        network[layer][node].setUB(maxUB);
      }

    }
  }

  double lengthSoFar = 1;
  double chisqSoFar = 0;
  double pvalue = 0;
  int pastSize = 0;


  // traverse the network in a breadth-first strategy
  for (int layer = ncols; layer >= 1; layer--) {
    for (size_t node = 0; node < network[layer].size(); node++) {
      pastSize = network[layer][node].getPastSize();

      for (int i = 0; i < pastSize; i++) {

        chisqSoFar = network[layer][node].getPastChisq(i);

        if (chisqSoFar + network[layer][node].getUB() < funchisq) {}
        else {

          lengthSoFar = network[layer][node].getPastLen(i);

          if (chisqSoFar + network[layer][node].getLB() >= funchisq) {
            pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
          }
          else {
            for (int k = 0; k < network[layer][node].getSize(); k++) {
              network[layer - 1][network[layer][node].getChildrenIndex(k)].addPastLen(
                  network[layer][node].getLengthToChildren(k) * lengthSoFar,
                  network[layer][node].getColChisqToChildren(k) + chisqSoFar
              );
            } // end child
          }
        }
      } // end pastChisq

    } //end node
  } // end layer

  pvalue /= marginal;

  return pvalue;
}


