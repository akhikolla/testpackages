//  EFT_DQP.h
//    Exact functional test implementaiton using dynamic and
//    quadratic programming
//
//  Created by Hien Nguyen on 7/24/2018.
//
//  Revision history:
//  2019-02-25 (Hien Nguyen): created the namespace DQP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.h of 2.4.7 to EFT_DQP.h to distinguish from other
//     implementations of the exact functional test.

#include "Node.h"
#include "define.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>

namespace DQP {
//convertToInt: convert Rs to equivalent integer
int convertToInt(vector<int> Rs);

// isMember: to find the node in a layer
// if the node exists in the layer, return its index
// otherwise, return -1

int isMember(int &eInt, vector<Node> &layer);

int searchHashTable(vector<vector< pair<int, int>>> &hashTable, int element) ;

// compute the weight between two nodes
double colChisq(vector<int> &Rs1, vector<int> &Rs2, int sum, const vector<int> & squares, double &COLMARGIN) ;

double colChisq(vector<int> Rs1, int &sum, const vector<int> & squares, double &COLMARGIN) ;

// compute the length from the current node to the end node
double length(vector<int> &Rs1, int &sum, int &layer, vector<int> &Cs, const vector<double> & factorials);

// compute the length between two nodes
double length(vector<int> &Rs1, vector<int> &Rs2, const vector<double> & factorials);

// compute the funchisq without the fixed marginals
double funchisqByCol(vector<vector<int>> &observedTable, vector<int> &CSUM, const vector<int> & squares, double &COLMARGIN);
    
// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void createNode(Node &node, vector<int> Rs, vector<int> &Cs, int layer, vector<int> &currRs, int &nrows, int sum1, int sum2,
                vector<int> &S, const int &i, const vector<int> & squares, const vector<double> & factorials, vector<Node> &Layer,
                double &COLMARGIN, vector<vector<pair<int, int>>> &hashTable);

double lower_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums, double &COLMARGIN);


double upper_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums, double &COLMARGIN);

// main program
double EFTNetwork(vector<vector<int>> observedTable);
}
