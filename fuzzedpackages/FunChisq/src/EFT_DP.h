//  EFT_DP.h
//    Exact functional test implementaiton using dynamic programming
//
//  Created by Hien Nguyen on 6/17/2018.
//
//  Revision history:
//  2019-02-25 (Hien Nguyen): created the namespace DP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.h of 2.4.6 to EFT_DP.h to distinguish from other
//     implementations of the exact functional test.

#include "Node.h"
#include "define.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>

namespace DP {
//convertToInt: convert Rs to equivalent integer
int convertToInt(vector<int> Rs);

// isMember: to find the node in a layer
// if the node exists in the layer, return its index
// otherwise, return -1

int isMember(int &eInt, vector<Node> &layer);

int searchHashTable(vector<vector< pair<int, int>>> &hashTable, int element) ;

// compute the weight between two nodes
double colChisq(Node &node, vector<int> &Rs2, int &sum, const vector<int> & squares, double &COLMARGIN);

// compute the length from the current node to the end node
double length(Node &node, int &sum, int &layer, vector<int> &Cs, const vector<double> & factorials);

// compute the length between two nodes
double length(Node &node, vector<int> &Rs2, const vector<double> & factorials);

// compute the funchisq without the fixed marginals
double funchisqByCol(vector<vector<int>> &observedTable, vector<int> &CSUM, const vector<int> & squares, double &COLMARGIN);

// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void createNode(Node &node, vector<int> &Cs, int &layer, vector<int> &currRs, int &nrows, int sum1, int sum2,
                vector<int> &S, const int &i, const vector<int> & squares, const vector<double> & factorials,
                vector<Node> &Layer, double &COLMARGIN, vector<vector<pair<int, int>>> &hashTable);

// main program
double EFTNetwork(vector<vector<int>> observedTable);

}
