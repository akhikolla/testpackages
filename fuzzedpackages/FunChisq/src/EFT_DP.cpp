//  EFT_DP.cpp
//    Exact functional test implementaiton using dynamic programming
//
//  Created by Hien Nguyen on 6/17/2018.
//
//  Revision history:
//  2019-02-25 (Hien Nguyen): insulate code into the namespace DP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.cpp of 2.4.6 to EFT_DP.cpp to distinguish from other
//     implementations of the exact functional test.
//

//#include "Node.h"
#include "EFT_DP.h"
#include "trimTable.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>
#include <algorithm>
#include <climits>

using namespace std;

//convertToInt: convert Rs to equivalent integer
int DP::convertToInt(vector<int> Rs) {
	std::sort(Rs.begin(), Rs.end());

	int eInt = 0;
	for (size_t x = 0; x < Rs.size(); x++) {
		eInt *= 127;
		eInt += Rs[x];
	}
	return eInt;
}


// isMember: to find the node in a layer
// if the node exists in the layer, return its index
// otherwise, return -1

int DP::isMember(int &eInt, vector<Node> &layer) {
	for (size_t x = 0; x < layer.size(); x++) {
		if (eInt == layer[x].getEquiInt()) return (int) x;
	}
	return -1;
}

int DP::searchHashTable(vector<vector< pair<int, int>>> &hashTable, int element) {
  int i = element % hashTable.size();
  for (size_t j = 0; j < hashTable[i].size(); j++) {
    if (hashTable[i][j].first == element) return hashTable[i][j].second;
  }
  return -1;
}

// compute the weight between two nodes
double DP::colChisq(Node &node, vector<int> &Rs2, int &sum, const vector<int> & squares, double &COLMARGIN) {
  if (sum  > 0) {
  	double colchisq = 0.0;
	  for (size_t x = 0; x < Rs2.size(); x++) {
  		colchisq += squares[node.getRsum().at(x) - Rs2[x]];
  	}
	  colchisq =  colchisq * COLMARGIN / sum;
	  return (colchisq);
  }
  else return 0;
}

// compute the length from the current node to the end node
double DP::length(Node &node, int &sum, int &layer, vector<int> &Cs, const vector<double> & factorials) {
	double length = factorials[sum];
	for (size_t x = 0; x < node.getRsum().size(); x++) {
		length /= factorials[node.getRsum().at(x)];
	}
	for (int x = 0; x < layer; x++) {
		length /= factorials[Cs[x]];
	}
	return (length);
}

// compute the length between two nodes
double DP::length(Node &node, vector<int> &Rs2, const vector<double> & factorials) {
	double length = 1.0;
	for (size_t x = 0; x < Rs2.size(); x++) {
		length /= factorials[node.getRsum().at(x)- Rs2[x]];
	}
	return (length);
}

// compute the funchisq without the fixed marginals
double DP::funchisqByCol(vector<vector<int>> &observedTable,
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

// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void DP::createNode(Node &node, vector<int> &Cs, int &layer, vector<int> &currRs, int &nrows, int sum1, int sum2,
                    vector<int> &S, const int &i, const vector<int> & squares, const vector<double> & factorials,
                    vector<Node> &Layer, double &COLMARGIN, vector<vector<pair<int, int>>> &hashTable) {
  if (i == nrows) {
    double len = DP::length(node, currRs, factorials);
    double colchisq = DP::colChisq(node, currRs, Cs[layer - 1], squares, COLMARGIN);

    int eInt = DP::convertToInt(currRs);
    int index = DP::searchHashTable(hashTable, eInt);

    // if the child node does not exist yet, insert it to the next layer as a new node
    if (index < 0) {
      Layer.push_back(Node(currRs, eInt));
      node.addChildLink((int) Layer.size() - 1, len, colchisq);
      //update hashTable by insertion
      hashTable[eInt % hashTable.size()].push_back(make_pair(eInt,Layer.size() - 1));
    }

    // if the child node already exists, add a new link from the current node to that child node
    else {
      node.addChildLink(index, len, colchisq);
    }

  }
  else {
    int lowerbound, upperbound;

    sum1 += (i > 0 ? node.getRsum()[i - 1] : 0);
    sum2 += (i > 0 ? currRs[i - 1]: 0);

    lowerbound = std::max(0, node.getRsum().at(i) - Cs[layer-1] + sum1 - sum2);
    upperbound = std::min(node.getRsum().at(i), (layer-2 >= 0 ? S[layer-2]: 0) - sum2);

    for (int x = lowerbound; x <= upperbound; x++) {
      currRs[i] = x;
      DP::createNode(node, Cs, layer, currRs, nrows, sum1, sum2, S, i + 1, squares, factorials, Layer, COLMARGIN, hashTable);
    }
  }
}

double DP::EFTNetwork(vector<vector<int>> observedTable) {

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

	double COLMARGIN = 1.0;
	for (int i = 0; i < ncols; i++) {
	  if (ColSums[i]  > 0 ) COLMARGIN *= ColSums[i];
	}


	vector<vector<Node>> network(ncols + 1);

	//create the start node
	network[ncols].push_back(Node(GlobalRowSums, 0));
	network[ncols][0].addPastLen(1.0, 0);


	// hash table
	int hashTableSize = 199;
	vector<vector< pair<int, int>>> hashTable(hashTableSize);

	// generate the network without Quadratic bounds
	vector<int> currRs(nrows);
	for (int layer = ncols; layer > 0; layer--) {
	  for (size_t n = 0; n < network[layer].size(); n++) {
	    DP::createNode(network[layer][n], ColSums, layer, currRs, nrows, 0, 0, S, 0, squares, factorials, network[layer - 1], COLMARGIN, hashTable);
	  }

	 // for (int x = 0; x < hashTableSize; x++) hashTable[x].clear(); //clear hashtable
	}


	// compute upperbound and lowerbound for Layer1
	for (size_t x = 0; x < network[1].size(); x++) {
		network[1][x].setLB(network[1][x].getColChisqToChildren(0));
		network[1][x].setUB(network[1][x].getColChisqToChildren(0));
		network[1][x].setLengthToEnd(network[1][x].getLengthToChildren(0));
	}

	// compute upperbound and lowerbound for higher layer
	double minLB = 0;
	double maxUB = 0;
	double tempBound, colchisq;

	for (int layer = 2; layer <= ncols; layer++) {
		for (size_t node = 0; node < network[layer].size(); node++) {

			network[layer][node].setLengthToEnd(DP::length(network[layer][node], S[layer - 1], layer, ColSums, factorials));

			minLB = INT_MAX; maxUB = 0;

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

	double lengthSoFar = 1;
	double chisqSoFar = 0;
	double pvalue = 0;
	int pastSize = 0;

	// compute the adjusted funchisq
	double funchisq = DP::funchisqByCol(observedTable, ColSums, squares,COLMARGIN);

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


