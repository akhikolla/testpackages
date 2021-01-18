
#include "Node.h"
#include "define.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>

using namespace std;

//convertToInt: convert Rs to equivalent integer
int convertToInt(vector<int> Rs);

// isMember: to find the node in a layer
// if the node exists in the layer, return its index
// otherwise, return -1

int isMember(int &eInt, vector<Node> &layer);

int searchHashTable(vector<vector< pair<int, int>>> &hashTable, int element) ;

// compute the weight between two nodes
double colChisq(vector<int> &Rs1, vector<int> &Rs2, int &sum, int squares[], int &COLMARGIN) ;
 
double colChisq(vector<int> &Rs1, int &sum, int squares[], int &COLMARGIN);

// compute the length from the current node to the end node
double length(vector<int> &Rs1, int &sum, int &layer, vector<int> &Cs, double factorials[]);

// compute the length between two nodes
double length(vector<int> &Rs1, vector<int> &Rs2, double factorials[]);

// compute the funchisq without the fixed marginals
double funchisqByCol(vector<vector<int>> &observedTable, vector<int> &CSUM, int squares[], int &COLMARGIN);

// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void createNode(Node &node,  vector<int> Rs, vector<int> &Cs, int layer, vector<int> &currRs, int nrows, int sum1, int sum2,
	vector<int> &S, int i, int squares[], double factorials[], vector<Node> &Layer, int &COLMARGIN,vector<vector<pair<int, int>>> &hashTable) ;

double lower_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums);
  
double upper_bound(int layer, vector <int> &Rsum, vector<int> &O_colsums);

// main program
double EFTNetwork(vector<vector<int>> observedTable);

