
#include "Node.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>
#include <algorithm>
#include <climits>

using namespace std;

//convertToInt: convert Rs to equivalent integer
int convertToInt(vector<int> Rs) {
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

int isMember(int &eInt, vector<Node> &layer) {
	for (size_t x = 0; x < layer.size(); x++) {
		if (eInt == layer[x].getEquiInt()) return (int) x;
	}
	return -1;
}

int searchHashTable(vector<vector< pair<int, int>>> &hashTable, int element) {
  int i = element % hashTable.size();
  for (size_t j = 0; j < hashTable[i].size(); j++) {
    if (hashTable[i][j].first == element) return hashTable[i][j].second;
  }
  return -1;
}

// compute the weight between two nodes
double colChisq(vector<int> &Rs1, vector<int> &Rs2, int &sum, int squares[], int &COLMARGIN) {
  if (sum  > 0) {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs2.size(); x++) {
      colchisq += squares[Rs1[x] - Rs2[x]];
    }
    colchisq =  colchisq * COLMARGIN / sum;
    colchisq = colchisq * Rs2.size() - COLMARGIN * sum;
    return (colchisq);
  }
  else return 0;
}

double colChisq(vector<int> Rs1, int &sum, int squares[], int &COLMARGIN) {
  if (sum  > 0) {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs1.size(); x++) {
      colchisq += squares[Rs1[x]];
    }
    colchisq =  colchisq * COLMARGIN / sum;
    colchisq = colchisq *Rs1.size() - COLMARGIN * sum;
    return (colchisq);
  }
  else return 0;
  
}
// compute the length from the current node to the end node
double length(vector<int> &Rs1, int &sum, int &layer, vector<int> &Cs, double factorials[]) {
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
double length(vector<int> &Rs1, vector<int> &Rs2, double factorials[]) {
  double length = 1.0;
  for (size_t x = 0; x < Rs2.size(); x++) {
    length /= factorials[Rs1[x]- Rs2[x]];
  }
  return (length);
}

// compute the funchisq without the fixed marginals
double funchisqByCol(vector<vector<int>> &observedTable, 
                       vector<int> &CSUM, int squares[], int &COLMARGIN) {
	double funchisq = 0.0;
	double colchisq = 0.0;
	for (size_t j = 0; j < observedTable[0].size(); j++) {
		colchisq = 0;
	 // if (CSUM[j] > 0) {
		  for (size_t i = 0; i < observedTable.size(); i++) {
			  colchisq += squares[observedTable[i][j]];
		  }
		  colchisq = colchisq * COLMARGIN / CSUM[j];
	 // }
		
	  funchisq += colchisq *observedTable.size() - COLMARGIN * CSUM[j];
	}
	return (funchisq);
}


double lower_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums, int &COLMARGIN)
{
  double lower_bound = 0;
  double temp;
  
  vector<int> U(Rsum);
  
  size_t nrows = Rsum.size();
  //size_t ncols = O_colsums.size();
  
  
  vector<size_t> order(nrows);
  for (size_t q = 0; q < nrows; ++q) {
    order[q] = q;
  }
  
  // sort U in increasing order
  sort(order.begin(), order.end(),
       [&U](size_t i1, size_t i2) {return U[i1] < U[i2]; });
  
  for (int l = layer - 1; l >= 0; l--) {
    // find lower bound for row l
    int runsum = 0;
    double e = O_colsums[l] / (double)nrows;
    for (size_t k = 0; k < nrows; ++k) {
      // accumulate the lower bound
      double xavg = (O_colsums[l] - runsum) / (double)(nrows - k);
      if (U[order[k]] < xavg) {
        temp = U[order[k]] - e;
        if (e>0) lower_bound += temp * temp *COLMARGIN / e;
        runsum += U[order[k]];
      }
      else {
        temp = xavg - e;
        if (e>0) lower_bound += (temp * temp * COLMARGIN/ e) * (nrows - k);
        break;
      }
    }
  }
  
  return lower_bound;
}


double upper_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums, int &COLMARGIN)
{
  double upper_bound = 0;
  double temp;
  vector<int> U(Rsum);
  
  //size_t ncols = O_colsums.size();
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
    int runsum = 0;
    double el = O_colsums[l] / (double)nrows;
    for (size_t k = 0; k < nrows; ++k) {
      // accumulate the lower bound
      int xmax = O_colsums[l] - runsum;
      if (U[order[k]] < xmax) {
        temp = U[order[k]] - el ;
        if (el>0) upper_bound += temp * temp * COLMARGIN / el;
        runsum += U[order[k]];
      }
      else if (xmax != 0) {
        temp = xmax - el;
        if (el>0) upper_bound += temp * temp * COLMARGIN/ el;
        runsum += xmax;
      }
      else {
        upper_bound += el * (nrows - k) * COLMARGIN;
        break;
      }
    }
  }
  
  return upper_bound;
  
}


// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void createNode(Node &node, vector<int> Rs, vector<int> &Cs, int layer, vector<int> &currRs, int &nrows, int sum1, int sum2,
                vector<int> &S, const int &i, int squares[], double factorials[], vector<Node> &Layer, int &COLMARGIN, vector<vector<pair<int, int>>> &hashTable) {
  if (i == nrows) {
    double len = length(Rs, currRs, factorials);
    int colchisq = colChisq(Rs, currRs, Cs[layer], squares, COLMARGIN);

    int eInt = convertToInt(currRs);
    int index = searchHashTable(hashTable, eInt);

    // if the child node does not exist yet, insert it to the next layer as a new node
    if (index < 0) {
      Layer.push_back(Node(currRs, eInt));
      node.addChildLink((int) Layer.size() - 1, len, colchisq);
      //update hashTable by insertion
      hashTable[eInt % hashTable.size()].push_back(make_pair(eInt,Layer.size() - 1));
      
      Layer[Layer.size() - 1].setMinPastChisq(node.getMinPastChisq() + colchisq);
      Layer[Layer.size() - 1].setMaxPastChisq(node.getMaxPastChisq() + colchisq);
      
      Layer[Layer.size() - 1].setLB(lower_bound(layer, currRs, Cs, COLMARGIN));
      Layer[Layer.size() - 1].setUB(upper_bound(layer, currRs, Cs, COLMARGIN));
      
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
      createNode(node, Rs, Cs, layer, currRs, nrows, sum1, sum2, S, i + 1, squares, factorials, Layer, COLMARGIN, hashTable);
    }
  }
}

double EFTNetwork(vector<vector<int>> observedTable) {

	int nrows = (int) observedTable.size();
	int ncols = (int) observedTable[0].size();

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
	
	for (int j = 0; j < ncols; j++) {
	  if (ColSums[j] == 0) {
	    for (int i = 0; i < nrows; i++) {
	      observedTable[i].erase(observedTable[i].begin() + j);
	    }
	    if (j < ncols-1) {
	      for (int k = j ; k < ncols -1; k++ ) ColSums[k] = ColSums[k + 1];
	    }
	    ncols--;
	  }
	}

	int squares[500];
	for (int x = 0; x < N; x++) squares[x] = x*x;

	double factorials[500];
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

	int COLMARGIN = 1;
	for (int i = 0; i < ncols; i++) {
	  if (ColSums[i]  > 0 ) COLMARGIN *= ColSums[i];
	}

	// compute the adjusted funchisq
	unsigned int funchisq = funchisqByCol(observedTable, ColSums, squares,COLMARGIN);

	vector<vector<Node>> network(ncols + 1);

	//create the start node
	network[ncols].push_back(Node(GlobalRowSums, 0));
	network[ncols][0].addPastLen(1.0, 0);
	network[ncols][0].setMaxPastChisq(0);
	network[ncols][0].setMinPastChisq(0);
	
	network[ncols][0].setLB(lower_bound(ncols, GlobalRowSums, ColSums, COLMARGIN));
	network[ncols][0].setUB(upper_bound(ncols, GlobalRowSums, ColSums, COLMARGIN));
	
	if (network[ncols][0].getLB() == funchisq) return 1;

	// hash table
	int hashTableSize = 199;
	vector<vector< pair<int, int>>> hashTable(hashTableSize);

	// generate the network
	vector<int> currRs(nrows);
	for (int layer = ncols; layer > 1; layer--) {
	  for (size_t n = 0; n < network[layer].size(); n++) {
	    if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() < funchisq) {
	      if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() >= funchisq) {
	        createNode(network[layer][n], network[layer][n].getRsum(), ColSums, layer-1, currRs, nrows, 0, 0, S, 0, squares, factorials, network[layer - 1], COLMARGIN, hashTable);
	      }
	    }
	  }

	 // for (int x = 0; x < hashTableSize; x++) hashTable[x].clear(); //clear hashtable
	}

	int colchisq;
	for (size_t x = 0; x < network[1].size(); x++) {
	  colchisq = colChisq(network[1][x].getRsum(), ColSums[0], squares, COLMARGIN);
	  network[1][x].setLB(colchisq);
	  network[1][x].setUB(colchisq);
	}

	double lengthSoFar = 1;
	int chisqSoFar = 0;
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
